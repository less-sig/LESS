#!/usr/bin/env python3
"""
Reads the NIST KAT files of `Utilities/KAT_Generation/KAT`, regenerates the
keys and signatures and verifies every signed message in them. This is
`Reference_Implementation/lib/nist/PQCgenKAT_sign.c`

usage:
    [LESS_PARAMS=cat3_220] python3 kat.py <path/to.rsp> [max_entries]
"""
import os
import sys

from config import LESS_CRYPTO_PUBLICKEYBYTES
from less import crypto_sign_keypair, crypto_sign, crypto_sign_open
from prng import init_randombytes


def parse(text: str) -> list[dict]:
    """ splits a `.rsp` file into its `count = ...` blocks """
    out = []
    for line in text.splitlines():
        key, _, value = line.partition("=")
        key, value = key.strip(), value.strip()
        if key == "count":
            out.append({})
        if out and key in ("count", "mlen", "smlen"):
            out[-1][key] = int(value)
        elif out and key in ("seed", "msg", "pk", "sk", "sm"):
            out[-1][key] = bytes.fromhex(value)
    return out


def check_entry(e: dict) -> str:
    init_randombytes(e["seed"])
    pk, sk = crypto_sign_keypair()
    if pk != e["pk"]:
        return "pk mismatch"
    if sk != e["sk"]:
        return "sk mismatch"

    sm = crypto_sign(e["msg"], sk)
    if sm != e["sm"] or len(sm) != e["smlen"]:
        return "sm mismatch"

    if crypto_sign_open(e["sm"], e["pk"]) != e["msg"]:
        return "crypto_sign_open rejected the signature"
    return ""


def main() -> int:
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        f"../../Utilities/KAT_Generation/KAT/PQCsignKAT_{LESS_CRYPTO_PUBLICKEYBYTES}.rsp")
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else None

    with open(path) as f:
        entries = parse(f.read())[:limit]

    print(path)
    print(f"public key bytes = {LESS_CRYPTO_PUBLICKEYBYTES}")
    print(f"entries = {len(entries)}")

    failed = 0
    for e in entries:
        msg = check_entry(e)
        print(f"count = {e['count']:<4} {msg or 'ok'}")
        failed += msg != ""

    print(f"verified = {len(entries) - failed}, failed = {failed}")
    return 0 if failed == 0 and entries else 4


if __name__ == "__main__":
    sys.exit(main())
