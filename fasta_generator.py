# padlock design pipeline for multiplexed assay with multiple probes per target in cDNA-based expression profiling
# Xiaoyan, 2017

from lib import checkinput

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
