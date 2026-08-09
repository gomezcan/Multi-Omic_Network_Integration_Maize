#!/bin/bash

### output ##
## IDR M TF deM Target
grep 'Zm' targets.2kb_TSS.*.txt | sed 's/targets.2kb_TSS.Q30.B73.//g' | sed $'s/.txt:/\t/g' | cut -f1,5 | sed $'s/_/\t/g' > $1
# version 
#targets.2kb_TSS.Q30.
#grep 'Zm' targets.2kb_TSS.*.txt | sed 's/targets.2kb_TSS.Q30.//g' | sed $'s/.txt:/\t/g' | cut -f1,5 | sed $'s/_/\t/g' > $1
