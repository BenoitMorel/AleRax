#/bin/sh

ask_yes_no() {
    while true; do
        read -p "$1 (y/N): " answer
        case $answer in
            [Yy]* ) return 0 ;;
            ""|[Nn]* ) return 1 ;;
            * ) echo "Please answer yes or no." ;;
        esac
    done
}

num_cores=$(nproc 2>/dev/null || echo 1)

# Example usage:
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release

if ask_yes_no "Build only AleRax binary (otherwise build everything)?"; then
    echo "Building AleRax..."
    make -j $num_cores alerax
else
    echo "Building everything..."
    make -j $num_cores
fi
