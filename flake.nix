{
  description = "Nuclear Measurement Utilities";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
    }:
    (flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config = {
            allowUnfree = true;
            cudaCapabilities = [ "12.0" ];
            cudaForwardCompat = false;
          };
        };

        version = "26.5.31";

        rootWithCuda = pkgs.root.overrideAttrs (old: {
          cmakeFlags = (old.cmakeFlags or [ ]) ++ [
            "-Dcuda=ON"
            "-DCMAKE_CUDA_ARCHITECTURES=120"
          ];
          buildInputs =
            (old.buildInputs or [ ])
            ++ (with pkgs.cudaPackages; [
              cuda_nvcc
              cuda_cudart
              cuda_cccl
            ]);
        });

        mkUtils =
          { cuda }:
          let
            root = if cuda then rootWithCuda else pkgs.root;
            cudaBuildInputs = pkgs.lib.optionals cuda (
              with pkgs.cudaPackages;
              [
                cuda_nvcc
                cuda_cudart
                cuda_cccl
              ]
            );
            cudaCFlag = pkgs.lib.optionalString cuda "-DAU_ROOFIT_BACKEND_CUDA=1";
          in
          pkgs.stdenv.mkDerivation {
            pname = "analysis-utilities" + pkgs.lib.optionalString cuda "-cuda";
            inherit version;

            src = ./.;

            nativeBuildInputs =
              with pkgs;
              [
                pkg-config
                cmake
              ]
              ++ pkgs.lib.optionals (!pkgs.stdenv.hostPlatform.isDarwin) [
                autoPatchelfHook
              ];

            buildInputs = [ root ] ++ cudaBuildInputs;

            cmakeFlags = pkgs.lib.optionals cuda [
              "-DAU_USE_CUDA=ON"
              "-DCMAKE_CUDA_ARCHITECTURES=120"
            ];

            postFixup = pkgs.lib.optionalString (!pkgs.stdenv.hostPlatform.isDarwin) ''
              for lib in $out/lib/*.so; do
                if [ -f "$lib" ]; then
                  patchelf --set-rpath "$out/lib:${root}/lib:${pkgs.stdenv.cc.cc.lib}/lib" "$lib" || true
                fi
              done
            '';

            setupHook = ./setup-hook.sh;
          };

        utils = mkUtils { cuda = false; };
        utilsCuda = mkUtils { cuda = true; };

        pythonPackage = pkgs.python3Packages.buildPythonPackage {
          pname = "analysis-utilities";
          inherit version;
          src = ./python;
          pyproject = true;

          nativeBuildInputs = [ pkgs.python3Packages.setuptools ];
          postPatch = ''
            substituteInPlace analysis_utilities/__init__.py \
              --replace-fail '@VERSION@' "${version}"
            substituteInPlace pyproject.toml \
              --replace-fail '@VERSION@' "${version}"
          '';
          propagatedBuildInputs = with pkgs.python3Packages; [
            numpy
            pandas
          ];
          doCheck = false;
        };

        # Drop-in .clangd for the gpu/ subdir so clangd can parse .cu files
        # against the same arch + version policy that nvcc will use.
        clangdConfigFile = (pkgs.formats.yaml { }).generate "dot-clangd" {
          CompileFlags.Add = [
            "--cuda-gpu-arch=sm_120"
            "--no-cuda-version-check"
          ];
          Diagnostics.Suppress = [
            "no_member"
            "nested_name_spec_non_tag"
            "typename_nested_not_found"
            "template_instantiate_undefined"
          ];
        };
      in
      {
        # CPU-batched build (default). Downstream flakes that just want the
        # speedup from doEval pull this.
        packages.default = utils;
        # CUDA build — RunFit() and friends call EvalBackend::Cuda(). Downstream
        # flakes that opt in pull utils.packages.${system}.cuda.
        packages.cuda = utilsCuda;
        packages.pythonPackage = pythonPackage;
        # Re-export the CUDA-overlaid ROOT so downstream consumers can pull it
        # without duplicating the overrideAttrs block. Plain pkgs.root is left
        # for downstream to grab directly from nixpkgs.
        packages.rootCuda = rootWithCuda;

        devShells.cpu = pkgs.mkShell {
          buildInputs = with pkgs; [
            root
            gnumake
            pkg-config
            clang-tools
            (python3.withPackages (
              python-pkgs: with python-pkgs; [
                numpy
                pandas
                pytest
                pythonPackage
              ]
            ))
          ];

          shellHook = ''
            export SHELL="${pkgs.bash}/bin/bash"
            echo "Development environment for working on the analysis utilities source (CPU)"
            export CPLUS_INCLUDE_PATH="$PWD/include''${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"
            export ROOT_INCLUDE_PATH="$PWD/include:${pkgs.root}/include"
            export LD_LIBRARY_PATH="$PWD/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
          '';
        };

        devShells.default = pkgs.mkShell {
          buildInputs =
            (with pkgs; [
              gnumake
              pkg-config
              clang-tools
              (python3.withPackages (
                python-pkgs: with python-pkgs; [
                  numpy
                  pandas
                  pytest
                  pythonPackage
                ]
              ))
            ])
            ++ [ rootWithCuda ]
            ++ (with pkgs.cudaPackages; [
              cuda_nvcc
              cuda_cudart
              cuda_cccl
            ]);

          shellHook = ''
            export SHELL="${pkgs.bash}/bin/bash"
            echo "Development environment for working on the analysis utilities source (CUDA, AU_ROOFIT_BACKEND_CUDA=1)"
            export NIX_CFLAGS_COMPILE="-DAU_ROOFIT_BACKEND_CUDA=1''${NIX_CFLAGS_COMPILE:+ $NIX_CFLAGS_COMPILE}"
            export CPLUS_INCLUDE_PATH="$PWD/include''${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"
            export ROOT_INCLUDE_PATH="$PWD/include:${rootWithCuda}/include"
            export LD_LIBRARY_PATH="$PWD/lib:/run/opengl-driver/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
            mkdir -p "$PWD/gpu"
            install -m 644 ${clangdConfigFile} "$PWD/gpu/.clangd"
          '';
        };
      }
    ))
    // {
      templates = {
        default = {
          path = ./templates/standard;
          description = "Standard analysis development environment.";
          welcomeText = ''
            Run `nix develop` to enter the development environment.
            If you have local libraries in include/src, use the included Makefile, and run your macros with root -l macro.cpp+.
          '';
        };
        standard = self.templates.default;
        python = {
          path = ./templates/python;
          description = "Python analysis development environment with analysis-utils and ML libraries.";
          welcomeText = ''
            Run `nix develop` to enter the development environment.
            The analysis-utils Python package and PlottingUtils bridge are available out of the box.
          '';
        };
      };
    };
}
