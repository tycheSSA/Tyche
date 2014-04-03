template<typename T>
void DiffusionWithTracking<T>::integrate(const double dt) {
  std::set<int> dirty_indices;
  const int n = get_species().size();
  for (int i = 0; i < n; ++i) {
    Species &s = *(get_species()[i]);
    const double step_length = calc_step_length(s, dt);
    const int n = s.mols.size();
    std::vector<int> copy_numbers(indices_to_consider.size(), 0);
    for (int j = 0; j < n; ++j) {
      const Vect3d r0 = s.mols.r[j];
      const Vect3d r = r0+step_length * Vect3d(norm(),norm(),norm());
      //const int old_cidx = s.grid->get_cell_index(r0);
      const int cidx = s.grid->get_cell_index(r);
      for (int l = 0; l < indices_to_consider.size(); l++) {
	if (cidx==indices_to_consider[l]) {
	  copy_numbers[l]++;
	  break;
	}
      }
      s.mols.r0[j] = r0;
      s.mols.r[j] = r;
    }
    for (int k = 0; k < indices_to_consider.size(); k++) {
      const int cidx = indices_to_consider[k];
      if (s.copy_numbers[cidx]!=copy_numbers[k]) {
	s.copy_numbers[cidx] = copy_numbers[k];
	dirty_indices.insert(cidx);
      }
    }
  }
  for (auto i : dirty_indices)
    nsm->recalc_priority(i);
}
