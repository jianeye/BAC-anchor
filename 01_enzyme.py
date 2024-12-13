import re
import sys
import os

def find_cutter_sites(genome_seq, cutter_seq="TTCGAA"):
    positions = [m.start() for m in re.finditer(cutter_seq, genome_seq)]
    distances = []
    if len(positions) > 1:
        distances = [positions[i+1] - positions[i] for i in range(len(positions)-1)]
    
    return positions, distances

if __name__ == "__main__":
    inp = sys.argv[1]
    oup = inp + ".no" 

    with open(inp, "r") as f:
        genome_seq = "".join([line.strip() for line in f.readlines() if not line.startswith(">")])

    positions, distances = find_cutter_sites(genome_seq)

    print(f"no: {len(positions)}")
    """
	if len(positions) > 1:
        print("dis (bp):")
        print(distances)
    else:
        print("na")
    """
    with open(oup, "w") as out_f:
        #out_f.write(f"no: {len(positions)}\n")
        if len(positions) > 1:
            #out_f.write("dis (bp):\n")
            out_f.write("\n".join(map(str, distances)) + "\n")
        else:
            out_f.write("na\n")
