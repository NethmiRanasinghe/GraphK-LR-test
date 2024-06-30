#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status

cpath=$(pwd)
overall_success=true

# Function to check command success
check_success() {
    if [ $? -ne 0 ]; then
        echo "Error: $1 failed"
        overall_success=false
    fi
}

# Create a directory for binaries
[ -d bin ] && rm -rf bin
mkdir bin

# Install and build wtdbg2
[ -d wtdbg2 ] && rm -rf wtdbg2
git clone https://github.com/ruanjue/wtdbg2
(cd wtdbg2 && make) || check_success "wtdbg2 build"
mv wtdbg2/kbm2 bin/ || check_success "moving kbm2"
rm -rf wtdbg2

# Install and build seqtk
[ -d seqtk ] && rm -rf seqtk
git clone https://github.com/lh3/seqtk.git
(cd seqtk && make) || check_success "seqtk build"
mv seqtk/seqtk bin/ || check_success "moving seqtk"
rm -rf seqtk

# Install and build seq2vec
[ -d seq2vec ] && rm -rf seq2vec
git clone https://github.com/anuradhawick/seq2vec.git
(cd seq2vec && sh build.sh) || check_success "seq2vec build"
mv seq2vec/build/seq2vec bin/ || check_success "moving seq2vec"
#rm -rf seq2vec

# Install and build seq2covvec
[ -d seq2covvec ] && rm -rf seq2covvec
cd bin
git clone https://github.com/anuradhawick/seq2covvec.git
(cd seq2covvec && sh build.sh) || check_success "seq2covvec build"
cd $cpath

case $1 in
  osx | macos)
    echo "OSX Build"
    # OSX Build
    clang++ splitreads.cpp -o graphk/oblr_utils/split -lz -std=c++11 -O3 || check_success "OSX splitreads build"
    ;;
  *)
    echo "Linux Build"
    # Linux Build
    g++ splitreads.cpp -o graphk/oblr_utils/split -lz -std=c++11 -O3 || check_success "Linux splitreads build"
    ;;
esac

if $overall_success; then
    echo "All installations and builds completed successfully!"
else
    echo "Some installations or builds failed. Please check the output for errors."
    exit 1
fi