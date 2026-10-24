"""Serve the local demonstrations over HTTP on this computer."""

import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()
    handler = partial(SimpleHTTPRequestHandler, directory=str(ROOT))
    with ThreadingHTTPServer(("127.0.0.1", args.port), handler) as server:
        print(f"Demo navigation: http://localhost:{args.port}/demos/index.html", flush=True)
        print(f"Mesh optimization flow: http://localhost:{args.port}/results/mentor_2026-10-23/optimization_flow/index.html", flush=True)
        print(f"Physical-time evolution: http://localhost:{args.port}/results/mentor_2026-10-23/p2_v2/index.html", flush=True)
        print("Press Ctrl+C to stop this preview.", flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    main()
