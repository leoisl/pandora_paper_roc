ids = list(range(4))
ranges = [(20000, 30000), (880000, 890000), (1850000, 1860000), (1380000, 1390000), (750000, 760000),
		  (4610000, 4620000), (380000, 390000), (4080000, 4090000)]
for id_1 in ids:
	for id_2 in ids:
		if id_1 < id_2:
			conditions = []
			for (start, end) in ranges:
				conditions.append(f"($2 >= {start} && $2 <= {end})")
			conditions_str = " || ".join(conditions)

			print(f"head -n 1 genome.{id_1}-SEP-genome.{id_2}.mummer.csv | cat - "
				  f"<(awk '{conditions_str}' genome.{id_1}-SEP-genome.{id_2}.mummer.csv) > genome.{id_1}-SEP-genome.{id_2}.mummer.filtered.csv")