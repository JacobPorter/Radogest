#!/bin/sh
mkdir virus
mkdir bacteria
mkdir eukaryotes
mkdir archaea
~/Applications/Radogest/utils/change_taxid.py 2 2157 2759 1> ./virus/1.taxid &
~/Applications/Radogest/utils/change_taxid.py 10239 2157 2759 1> ./bacteria/1.taxid &
~/Applications/Radogest/utils/change_taxid.py 10239 2157 2 1> ./eukaryotes/1.taxid &
~/Applications/Radogest/utils/change_taxid.py 2 10239 2759 1> ./archaea/1.taxid &
