ids = list(range(4))
ranges = [(0, 50000), (2370000, 2720000)]
for id_1 in ids:
	for id_2 in ids:
		if id_1 < id_2:
			conditions = []
			for (start, end) in ranges:
				conditions.append(f"($2 >= {start} && $2 <= {end})")
			conditions_str = " || ".join(conditions)

			print(f"head -n 1 genome.{id_1}-SEP-genome.{id_2}.mummer.csv | cat - "
				  f"<(awk '{conditions_str}' genome.{id_1}-SEP-genome.{id_2}.mummer.csv) > genome.{id_1}-SEP-genome.{id_2}.mummer.filtered.csv")