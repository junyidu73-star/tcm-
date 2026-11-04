from types import SimpleNamespace

from research_reliability import _install_database_retry_patch


def test_database_retry_retries_transient_failures(monkeypatch):
    monkeypatch.setattr("research_reliability.time.sleep", lambda _: None)

    class FakeEngine:
        def __init__(self):
            self.attempts = 0
            self.job = {}
            self.config = {}

        def get(self, url, payload=None):
            self.attempts += 1
            if self.attempts < 3:
                raise TimeoutError("upstream timeout")
            return b"{}", "application/json"

    module = SimpleNamespace(ResearchEngine=FakeEngine)
    _install_database_retry_patch(module)
    engine = module.ResearchEngine()
    assert engine.get("https://example.org/query") == (b"{}", "application/json")
    assert engine.attempts == 3
