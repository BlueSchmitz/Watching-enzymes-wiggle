# usage: python scale_all_residues.py > processed_scaled_all.top
# adds "_" to all residue names in processed.top

section = None

with open("../processed.top") as f:
    for line in f:
        stripped = line.strip()
        # Track section
        if stripped.startswith("["):
            section = stripped.lower()
            print(line, end="")
            continue
        if stripped == "" or stripped.startswith(";"):
            print(line, end="")
            continue
        parts = line.split()
        # Only process atoms section
        if section == "[ atoms ]":
            try:
                resnum = int(parts[2])
            except:
                print(line, end="")
                continue
            parts[1] = parts[1] + "_"
        print("  ".join(parts))
