let
  pkgs = (import <nixpkgs> {}).pkgsCross.armv7l-hf-multiplatform;
in
{ 
    # pkgs ? import <nixpkgs> {
    #     localSystem = {
    #     gcc.arch = "znver4";
    #     gcc.tune = "znver4";
    #     system = "x86_64-linux";
    #   };
    # } 
}:
pkgs.gccStdenv.mkDerivation {
  name = "LESS";
  src = ./.;
  nativeBuildInputs = with pkgs; [ 
    file 
    pkg-config
    gnumake 
    cmake 
    gcc
    gdb
    glibc
    openssl
    gcc 
    clang
  ];

  buildInputs = with pkgs; [
    zlib
  ];
}#){}
