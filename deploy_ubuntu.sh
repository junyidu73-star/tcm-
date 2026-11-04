#!/usr/bin/env bash
set -euo pipefail

DOMAIN="${1:-}"
ADMIN_EMAIL="${2:-}"
SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP_ROOT="/opt/zhiheng"

fail() { printf '错误：%s\n' "$1" >&2; exit 1; }

[[ "${EUID}" -eq 0 ]] || fail "请使用 sudo 运行：sudo bash deploy_ubuntu.sh 域名 管理员邮箱"
[[ -n "$DOMAIN" && -n "$ADMIN_EMAIL" ]] || fail "用法：sudo bash deploy_ubuntu.sh qklkwznd.colocrossing.cloud you@example.com"
[[ "$DOMAIN" =~ ^[A-Za-z0-9.-]+$ ]] || fail "域名格式不正确。"
[[ "$ADMIN_EMAIL" =~ ^[^[:space:]@]+@[^[:space:]@]+\.[^[:space:]@]+$ ]] || fail "管理员邮箱格式不正确。"
[[ -f "$SOURCE_DIR/server.py" && -f "$SOURCE_DIR/requirements.txt" ]] || fail "请从解压后的平台根目录运行脚本。"
[[ ! -e "$APP_ROOT" ]] || fail "$APP_ROOT 已存在。为避免覆盖现有部署，脚本停止；请先检查服务器上的目录。"

source /etc/os-release
[[ "${ID:-}" == "ubuntu" && "${VERSION_ID:-}" == "22.04" ]] || fail "此安装脚本针对 Ubuntu 22.04。"
command -v python3.13 >/dev/null 2>&1 || fail "当前研究引擎依赖 CPython 3.13 字节码。请先在完成服务器环境审查后安装并确认 python3.13；本脚本不会擅自添加第三方软件源。"
getent ahostsv4 "$DOMAIN" >/dev/null || fail "$DOMAIN 当前无法解析到 IPv4。请先在 VPS 控制面板核对 DNS。"

printf '安装系统依赖……\n'
apt-get update
DEBIAN_FRONTEND=noninteractive apt-get install -y python3 python3-venv python3-pip nginx certbot rsync openssl
if command -v ufw >/dev/null 2>&1 && ufw status | grep -q 'Status: active'; then
    ufw allow 80/tcp
    ufw allow 443/tcp
fi

install -d -m 0755 /var/www/letsencrypt/.well-known/acme-challenge
cat > /etc/nginx/sites-available/zhiheng-platform <<EOF
server {
    listen 80;
    server_name ${DOMAIN};
    location ^~ /.well-known/acme-challenge/ {
        root /var/www/letsencrypt;
    }
    location / { return 404; }
}
EOF
ln -sfn /etc/nginx/sites-available/zhiheng-platform /etc/nginx/sites-enabled/zhiheng-platform
nginx -t
systemctl reload nginx || systemctl start nginx

printf '申请 HTTPS 证书……\n'
certbot certonly --non-interactive --agree-tos --webroot \
    --webroot-path /var/www/letsencrypt --email "$ADMIN_EMAIL" -d "$DOMAIN"

printf '安装平台文件和 Python 依赖……\n'
if ! id -u zhiheng >/dev/null 2>&1; then
    useradd --system --home-dir "$APP_ROOT" --shell /usr/sbin/nologin zhiheng
fi
install -d -o zhiheng -g zhiheng -m 0750 "$APP_ROOT"
rsync -a \
    --exclude='.venv/' --exclude='users.sqlite3' --exclude='*.sqlite3' \
    --exclude='__pycache__/' --exclude='runs/' --exclude='*.log' --exclude='*.bak' \
    --exclude='pytest-cache-files-*/' --exclude='tests/tmp*/' \
    "$SOURCE_DIR/" "$APP_ROOT/"
"$(command -v python3.13)" -m venv "$APP_ROOT/.venv"
"$APP_ROOT/.venv/bin/python" -m pip install --no-cache-dir -r "$APP_ROOT/requirements.txt"
chown -R zhiheng:zhiheng "$APP_ROOT"

INVITE_CODE="$(openssl rand -hex 16)"
install -d -m 0750 /etc/zhiheng
umask 077
cat > /etc/zhiheng/platform.env <<EOF
ZHIHENG_PORT=8765
ZHIHENG_HTTPS=1
ZHIHENG_INVITE_CODE=${INVITE_CODE}
ZHIHENG_TRUSTED_PROXY_IPS=127.0.0.1,::1
EOF
printf '域名：https://%s\n邀请码：%s\n' "$DOMAIN" "$INVITE_CODE" > /root/zhiheng-invite-code.txt
chmod 0600 /etc/zhiheng/platform.env /root/zhiheng-invite-code.txt

cat > /etc/systemd/system/zhiheng-platform.service <<'EOF'
[Unit]
Description=Zhiheng Medical Mechanism Research Platform
After=network-online.target
Wants=network-online.target

[Service]
Type=simple
User=zhiheng
Group=zhiheng
WorkingDirectory=/opt/zhiheng
EnvironmentFile=/etc/zhiheng/platform.env
ExecStart=/opt/zhiheng/.venv/bin/python /opt/zhiheng/server.py
Restart=on-failure
RestartSec=5
NoNewPrivileges=true
PrivateTmp=true
ProtectSystem=full
ProtectHome=true

[Install]
WantedBy=multi-user.target
EOF

cat > /etc/nginx/conf.d/zhiheng-rate-limit.conf <<'EOF'
limit_req_zone $binary_remote_addr zone=zhiheng_auth:10m rate=12r/m;
EOF
cat > /etc/nginx/snippets/zhiheng-proxy.conf <<'EOF'
proxy_pass http://127.0.0.1:8765;
proxy_http_version 1.1;
proxy_set_header Host $host;
proxy_set_header X-Real-IP $remote_addr;
proxy_set_header X-Forwarded-For $remote_addr;
proxy_set_header X-Forwarded-Proto https;
proxy_read_timeout 300s;
proxy_send_timeout 300s;
EOF
cat > /etc/nginx/sites-available/zhiheng-platform <<EOF
server {
    listen 80;
    server_name ${DOMAIN};
    location ^~ /.well-known/acme-challenge/ {
        root /var/www/letsencrypt;
    }
    location / { return 301 https://\$host\$request_uri; }
}

server {
    listen 443 ssl;
    server_name ${DOMAIN};
    ssl_certificate /etc/letsencrypt/live/${DOMAIN}/fullchain.pem;
    ssl_certificate_key /etc/letsencrypt/live/${DOMAIN}/privkey.pem;
    ssl_protocols TLSv1.2 TLSv1.3;
    client_max_body_size 2m;

    location ~ ^/api/auth/(login|register)$ {
        limit_req zone=zhiheng_auth burst=5 nodelay;
        include /etc/nginx/snippets/zhiheng-proxy.conf;
    }
    location / {
        include /etc/nginx/snippets/zhiheng-proxy.conf;
    }
}
EOF
cat > /etc/letsencrypt/renewal-hooks/deploy/reload-zhiheng-nginx <<'EOF'
#!/bin/sh
systemctl reload nginx
EOF
chmod 0755 /etc/letsencrypt/renewal-hooks/deploy/reload-zhiheng-nginx

systemctl daemon-reload
systemctl enable --now zhiheng-platform
nginx -t
systemctl reload nginx
if systemctl list-unit-files certbot.timer --no-legend 2>/dev/null | grep -q certbot.timer; then
    systemctl enable --now certbot.timer
fi

printf '\n部署完成：https://%s\n' "$DOMAIN"
printf '邀请码已写入 /root/zhiheng-invite-code.txt（仅 root 可读）；请私下发给老师。\n'
printf '同时在 ColoCrossing 防火墙/安全组开放 TCP 80 和 443；不要开放 8765。\n'
printf '服务状态：systemctl status zhiheng-platform\n日志：journalctl -u zhiheng-platform -f\n'
