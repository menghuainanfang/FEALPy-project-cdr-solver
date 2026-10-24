"""Check offline page links, embedded JSON and inline JavaScript syntax."""

import argparse
from html.parser import HTMLParser
import json
from pathlib import Path
import subprocess
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
PAGES = ["index.html", "demos/index.html", "results/e_variational_mesh.html",
         "results/f_poisson_comparison.html",
         "results/mentor_2026-10-23/e_convergence/e_variational_mesh.html",
         "results/mentor_2026-10-23/burgers/h_burgers.html"]


class Page(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links, self.scripts, self.current = [], [], None

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        for key in ("href", "src"):
            if key in attrs:
                self.links.append(attrs[key])
        if tag == "script":
            self.current = {"type": attrs.get("type", "javascript"), "text": ""}

    def handle_data(self, data):
        if self.current is not None:
            self.current["text"] += data

    def handle_endtag(self, tag):
        if tag == "script" and self.current is not None:
            self.scripts.append(self.current)
            self.current = None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--node", default="node")
    args = parser.parse_args()
    rows = []
    for relative in PAGES:
        path = ROOT / relative
        page = Page()
        page.feed(path.read_text(encoding="utf-8"))
        local = 0
        for link in page.links:
            target = urlsplit(link)
            if target.scheme or not target.path:
                continue
            if not (path.parent / unquote(target.path)).resolve().exists():
                raise AssertionError(f"Broken local link: {relative} -> {link}")
            local += 1
        for script in page.scripts:
            if script["type"] == "application/json":
                json.loads(script["text"])
            else:
                subprocess.run([args.node, "--check"], input=script["text"], text=True,
                               encoding="utf-8", capture_output=True, check=True)
        rows.append({"page": relative, "local_links": local, "inline_scripts": len(page.scripts),
                     "status": "static_checks_passed"})
    print(json.dumps({"pages": rows, "browser_interactions": "not_verified"}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
