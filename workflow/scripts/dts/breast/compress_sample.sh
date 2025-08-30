#!/bin/bash

FILE_PATH=$1
bgzip -f -c "$FILE_PATH" > "$FILE_PATH.gz"
tabix -p bed "$FILE_PATH.gz"
