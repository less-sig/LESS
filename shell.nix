with import <nixpkgs> {};
let
  my-python = pkgs.python3;
  python-with-my-packages = my-python.withPackages (p: with p; [
    scipy
    sage
	python-lsp-server
    pythonPackages.pandas
    pythonPackages.scipy
    pythonPackages.jupyter
    pythonPackages.pycrypto
  ]);
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
gccStdenv.mkDerivation {
  name = "LESS";
  src = ./.;
  buildInputs = with pkgs; [ 
    python-with-my-packages
    gnumake 
    cmake 
    clang
    clang-tools
    gcc
    gdb
    valgrind
    massif-visualizer
    openssl
    libxslt
    linuxKernel.packages.linux_6_6.perf
  ];
}
