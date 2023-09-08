# padlock design pipeline for multiplexed assay with multiple probes per target in cDNA-based expression profiling
# Xiaoyan, 2017

from lib import checkinput
#from lib import screenseq
#from lib import formatrefseq
#from lib import parblast
#from lib import readblast
#from lib import createoutput
#from lib import distributeprobes
#from lib import finalizeprobes


if __name__ == "__main__":
    try:
        # get keyboard inputs and prepare sequences
        designpars, outpars, genepars, designinput = \
            checkinput.getdesigninput()

        print("All finished!")

    except:
        import sys
        print (sys.exc_info()[0])
        import traceback
        print (traceback.format_exc())
    finally:
        print("Press Enter to continue ...")
        input()
