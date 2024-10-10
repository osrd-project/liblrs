import sys
import tomllib
import json

rust = tomllib.load(open("Cargo.toml", "rb"))
python = tomllib.load(open("python/Cargo.toml", "rb"))
wasm = json.load(open("wasm/package.json", "rb"))

rust_version = rust["package"]["version"]
python_version = python["package"]["version"]
wasm_version = wasm["version"]
if rust_version == python_version == wasm_version:
    print("all versions are the same")
    sys.exit(0)
else:
    print(f"version mismatch\n"
          f"  rust: {rust_version}\n"
          f"  python: {python_version}\n"
          f"  wasm: {wasm_version}\n")

    sys.exit(1)
