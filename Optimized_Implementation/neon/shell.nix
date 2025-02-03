let
  my-python = pkgs.python3;
  python-with-my-packages = my-python.withPackages (p: with p; [
    scipy
    sage
	python-lsp-server
    pythonPackages.pandas
    pythonPackages.numpy
    pythonPackages.scipy
    pythonPackages.jupyter
    pythonPackages.pycrypto
    pythonPackages.galois
    pythonPackages.bitarray
    pythonPackages.pycryptodome 
  ]);
  pkgs = (import <nixpkgs> {}).pkgsCross.aarch64-multiplatform;
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
#pkgs.pkgsStatic.callPackage ({ mkShell, zlib, pkg-config, file }: mkShell {
  name = "LESS";
  src = ./.;
  nativeBuildInputs = with pkgs; [ 
    file 
    pkg-config
    #python-with-my-packages
    gnumake 
    cmake 
    #clang
    #clang-tools
    #llvm
    gcc
    gdb
    glibc
    #valgrind
    #massif-visualizer
    openssl
    #libxslt
    #linuxKernel.packages.linux_6_6.perf
  ];

  buildInputs = with pkgs; [
    zlib
  ];
}#){}
