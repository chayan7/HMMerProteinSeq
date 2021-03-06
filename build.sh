#!/bin/sh

unset DYLD_LIBRARY_PATH
if ! echo "$0" | grep '\.sh$' > /dev/null; then
    printf 'Please run using "bash" or "sh", but not "." or "source"\\n' >&2
    return 1
fi


conv="$(conda -V 2>&1)"
printf "$conv"
condaExec="$(echo "$conv" | grep 'not' >&1)"


if [ "$condaExec" != "" ];then
	printf "\\nPlease install Conda properly and run this again, Thanks. \\n"
  exit 2
fi
if [ "$condaExec" = "" ];then
	printf "\\nAdding HMMerProteinSeq environment with all dependencies. \\n"
fi



THIS_DIR=$(DIRNAME=$(dirname "$0"); cd "$DIRNAME"; pwd)
THIS_FILE=$(basename "$0")
THIS_PATH="$THIS_DIR/$THIS_FILE"
ENV_PATH="$THIS_DIR/env"
MAC_FILE="$ENV_PATH/eHMM.txt"
LIN_FILE="$ENV_PATH/elinHMM.txt"

if [ "$(uname)" = "Linux" ];then
	conda create --name eHMM -f --file "$LIN_FILE"
	printf "\\nPlease type the following command from any terminal to activate FlaGs environment and run:\\n"
	printf ">> conda activate eHMM\\n"
	printf "\\nTo deactivate please type:\\n"
	printf ">> conda deactivate\\n"
	printf "\\n"
fi
if [ "$(uname)" = "Darwin" ]; then
	conda create --name eHMM -f --file "$MAC_FILE"
	printf "\\nPlease type the following command from any terminal to activate FlaGs environment and run:\\n"
	printf ">> conda activate eHMM\\n"
	printf "\\nTo deactivate please type:\\n"
	printf ">> conda deactivate\\n"
	printf "\\n"
fi

printf "One can run HMMerProteinSeq.py by activating eHMM environment from any terminal.\\n"
printf "\\nInstallation complete, Thanks. \\n"
exit 2
