"""Local dev server for the app: like `python -m http.server 8000`, but tells the browser never to cache, so an
edited script or data/reaction_rules.js is always the one that runs (no stale pages after a change).
    python tools/serve.py          (then open http://localhost:8000)
"""
import http.server, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 8000

class NoCache(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *a, **k):
        super().__init__(*a, directory=ROOT, **k)
    def end_headers(self):
        self.send_header("Cache-Control", "no-store, must-revalidate")
        self.send_header("Expires", "0")
        super().end_headers()

if __name__ == "__main__":
    print(f"serving {ROOT} at http://localhost:{PORT} (no caching)")
    http.server.ThreadingHTTPServer(("", PORT), NoCache).serve_forever()
