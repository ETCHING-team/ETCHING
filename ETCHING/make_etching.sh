echo "#!/usr/bin/env bash"
echo

version=$(cat ../version)
echo "VERSION=\"$version\""
echo
cat pipeline.sh
