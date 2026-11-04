# 部署到 ColoCrossing Ubuntu 22.04

本部署包使用 Nginx 提供 HTTPS；Python 应用继续绑定 `127.0.0.1:8765`，不会直接暴露应用端口。安装脚本会配置 systemd、HTTPS 证书、邀请码注册和反向代理客户端 IP 识别。**启动入口为 `server.py`，不是 Streamlit 主文件。**

当前研究引擎依赖 CPython 3.13 编译字节码文件 `research_engine.recovered.pyc`。部署脚本要求系统中已存在 `python3.13`，否则会停止；它不会为安装 Python 而更改现有软件源。请先完成服务器只读环境审查，再单独决定如何提供 Python 3.13。

## 部署前

1. 在 VPS 控制面板确认实例公网 IPv4。
2. 确认 `qklkwznd.colocrossing.cloud` 解析到该公网 IP。若不能使用此主机名，先用自己控制的域名添加 A 记录。
3. 在 ColoCrossing 防火墙/安全组开放 TCP 80、443 和实际 SSH 端口。不要开放 TCP 8765。
4. 准备一个可接收证书到期提醒的邮箱。

## 安装

将 GitHub 仓库克隆到 VPS，或将完整项目包上传到 VPS。确认项目根目录内存在 `server.py`、`requirements.txt` 和 `research_engine.recovered.pyc`，然后：

```bash
cd 知衡医学机制发现平台_最终版
sudo bash deploy_ubuntu.sh qklkwznd.colocrossing.cloud 你的邮箱@example.com
```

安装器会检查 Ubuntu 22.04、域名解析和必需文件；安装 Python/Nginx/Certbot 依赖、创建 `/opt/zhiheng` 服务目录，申请 Let's Encrypt 证书，并创建 systemd 服务。申请证书时公网 TCP 80 必须能到达 VPS。证书续期后 Nginx 会自动重载。

如果服务器已启用 UFW，安装器会放行 80/443；若 UFW 未启用，请在确认 SSH 规则正确后再自行启用。ColoCrossing 控制面板的网络防火墙也必须放行 80/443。

## 邀请老师

安装完成后访问 `https://qklkwznd.colocrossing.cloud`。首次访问选择创建账号，填写图形验证码和邀请码。邀请码保存在服务器 `/root/zhiheng-invite-code.txt`，由管理员通过私密渠道交给老师。邀请码不会嵌入网页代码或 Git 仓库。

如要停止新注册，在 `/etc/zhiheng/platform.env` 添加 `ZHIHENG_REGISTRATION_DISABLED=1` 并执行 `sudo systemctl restart zhiheng-platform`。重启会使所有人重新登录，并清除内存中的模型 API 配置；用户需再次设置 API Key。

常用运维命令：

```bash
sudo systemctl status zhiheng-platform
sudo journalctl -u zhiheng-platform -f
sudo nginx -t
sudo systemctl reload nginx
sudo ufw status
```

## 数据与安全

- 用户密码哈希写入 `/opt/zhiheng/users.sqlite3`；请定期备份此文件及 `runs/`。
- 模型 API Key 和登录会话保存在应用进程内存，服务重启后需要重新登录并填写 API Key。
- `ZHIHENG_HTTPS=1` 启用 Secure 登录 Cookie。Nginx 只代理本机应用，并把客户端 IP 通过受信任的本机代理头传给应用；登录与注册另有 Nginx 按真实客户端 IP 的请求速率限制。
- 邀请码允许持有者创建账号。只通过私密渠道发送；需要作废时更改环境文件中的邀请码并重启服务。
- 本版适合小范围课题组试用，不是经过第三方安全审计的医疗服务。网页结论属于科研候选，不构成诊疗建议。
