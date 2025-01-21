let
  pkgs = import <nixpkgs> {
    localSystem = "x86_64-linux";
    crossSystem = "riscv64-linux";
  };
in
pkgs.gccStdenv.mkDerivation {
  name = "LESS-riscv";
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
  ];

  buildInputs = with pkgs; [
    zlib
  ];
}#){}
