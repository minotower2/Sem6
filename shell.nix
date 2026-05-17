{ pkgs ? import <nixpkgs> {} }: 
    pkgs.mkShell {
    
    nativeBuildInputs = [
      pkgs.qt5.wrapQtAppsHook
      pkgs.makeWrapper
      pkgs.openssl
      pkgs.qtcreator
    ];

    shellHook = ''
      setQtEnvironment=$(mktemp)
      random=$(openssl rand -base64 20 | sed "s/[^a-zA-Z0-9]//g")
      makeWrapper "$(type -p sh)" "$setQtEnvironment" "''${qtWrapperArgs[@]}" --argv0 "$random"
      sed "/$random/d" -i "$setQtEnvironment"
      source "$setQtEnvironment"
    '';
  }
