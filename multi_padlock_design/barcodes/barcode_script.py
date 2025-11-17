import os
import time

from multi_padlock_design.barcodes import barcode as bar


def main():
    start_time = time.time()

    savedir = "13-mer-constrained_panel"

    os.makedirs(savedir, exist_ok=True)

    candidates = bar.make_constrained_library(13)

    #greedy_panel = bar.multi_restart_greedy(n = 10, d = 4, seed = 50, candidates = candidates)
    greedy_panel, distribution = bar.parallel_greedy(
        n = 13,
        d = 5,
        restarts = 20,
        candidates = candidates
    )
    #with open("/nemo/lab/znamenskiyp/home/users/colasa/code/multi_padlock_design/constrained_panel/1495_d4_10bp_barcodes.txt", "r") as f:
    #    greedy_panel = [line.strip().split(",")[1] for line in f if line.strip()]


    greedytime = time.time()
    print(f'elapsed for greedy: {start_time-greedytime:.2f} seconds')
    print(distribution)


    with open(os.path.join(savedir, f"{len(greedy_panel)}_d4_10bp_barcodes.txt"), "w") as f:
        for i, code in enumerate(greedy_panel):
            f.write(f"{i},{code}\n")


    #    improved_greedy, distribution = bar.parallel_local_improvement(greedy_panel, candidates, d = 4, max_passes = 5, processes = 8, base_seed=123 )

        #    print(distribution)

    #    improtime = time.time()
        #    print(f'elapsed for improved: {greedytime-improtime:.2f} seconds')

    #   with open(os.path.join(savedir, f"{len(improved_greedy)}_d4_10bp_barcodes_improved.txt"), "w") as f:
    #        for i, code in enumerate(improved_greedy):
    #            f.write(f"{i},{code}\n")


if __name__=="__main__":
    main()