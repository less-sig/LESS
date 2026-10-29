with import <nixpkgs> {};
{ pkgs ? import <nixpkgs> {} }:

stdenv.mkDerivation {
  name = "less_python";
  src = ./.;

  buildInputs = with pkgs; [
    python3
    python3Packages.sympy
    python3Packages.numpy
    python3Packages.setuptools
    python3Packages.virtualenv
    python3Packages.pylint
    python3Packages.matplotlib
    python3Packages.pycryptodome
  ];

  shellHook = ''
    SOURCE_DATE_EPOCH=$(date +%s)
    VENV=.venv
    if test ! -d $VENV; then
      virtualenv $VENV
    fi
    source ./$VENV/bin/activate
  '';
}
