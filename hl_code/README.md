This is the code for hub labeling, /chl/ is the folder containing our implemention of SHL. /sspexpcodes/ contains all other 3 algorithms.

/test is a shell script for the experiment. It tries to generate the weighted graph if it is not in folder /testdata/. To run the script to test mondial, the command is './test ./data/mondial.txt'



Note that in our SHL, we print Avg_label directly, while other indexes'print Avg_label-1​, so we add 1 to their result.

