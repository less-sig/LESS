with import <nixpkgs> {};
let
  my-python = pkgs.python3;
  python-with-my-packages = my-python.withPackages (p: with p; [
    scipy
    sage
	python-lsp-server
    pandas
    numpy
    scipy
    jupyter
    pycrypto
    galois
    bitarray
    pycryptodome 
  ]);
in
  { 
  }:
gccStdenv.mkDerivation {
  name = "LESS";
  src = ./.;
  buildInputs = with pkgs; [ 
    python-with-my-packages
    gnumake 
    cmake 
    clang
    clang-tools
    llvm
    gcc
    pkg-config
    openssl
    cargo-nextest
    cargo 
    rustc
  ] ++ (lib.optionals pkgs.stdenv.isLinux ([
    gdb
    valgrind
    massif-visualizer
      libxslt
      perf
      #linuxKernel.packages.linux_6_6.perf
  ]));
}
