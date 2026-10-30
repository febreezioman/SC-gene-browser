#!/usr/bin/env python3
"""
Static file server with HTTP Range support (igv.js needs it for BAM/FASTA).
`python -m http.server` ignores Range, so it cannot serve indexed files.

Serves the project root: open http://localhost:8000/web/

Usage: python3 scripts/serve.py [--port 8000] [--bind 127.0.0.1]
"""
from __future__ import annotations

import argparse
import os
import re
from functools import partial
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RANGE_RE = re.compile(r"bytes=(\d*)-(\d*)$")


class RangeHandler(SimpleHTTPRequestHandler):
    extensions_map = {**SimpleHTTPRequestHandler.extensions_map,
                      ".bam": "application/octet-stream", ".bai": "application/octet-stream",
                      ".fai": "text/plain", ".gtf": "text/plain", ".fasta": "text/plain"}

    def end_headers(self):
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("Cache-Control", "no-cache")
        super().end_headers()

    def do_GET(self):
        m = RANGE_RE.match(self.headers.get("Range", ""))
        path = self.translate_path(self.path)
        if not m or not os.path.isfile(path):
            return super().do_GET()
        size = os.path.getsize(path)
        start_s, end_s = m.groups()
        if start_s:
            start, end = int(start_s), min(int(end_s) if end_s else size - 1, size - 1)
        else:  # suffix range: last N bytes
            start, end = max(size - int(end_s), 0), size - 1
        if start >= size or start > end:
            self.send_response(HTTPStatus.REQUESTED_RANGE_NOT_SATISFIABLE)
            self.send_header("Content-Range", f"bytes */{size}")
            self.end_headers()
            return
        length = end - start + 1
        self.log_message("range %s %d-%d (%.1f MB)", os.path.basename(path), start, end, length / 1e6)
        self.send_response(HTTPStatus.PARTIAL_CONTENT)
        self.send_header("Content-Type", self.guess_type(path))
        self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
        self.send_header("Content-Length", str(length))
        self.end_headers()
        with open(path, "rb") as f:
            f.seek(start)
            while length > 0:
                chunk = f.read(min(1 << 20, length))
                if not chunk:
                    break
                self.wfile.write(chunk)
                length -= len(chunk)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--bind", default="127.0.0.1")
    args = ap.parse_args()
    server = ThreadingHTTPServer((args.bind, args.port), partial(RangeHandler, directory=str(ROOT)))
    print(f"Serving {ROOT} at http://{args.bind}:{args.port}/web/", flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()
