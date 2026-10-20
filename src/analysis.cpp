double calc_separation(const std::array<double, 3>& pos_i, const std::array<double, 3>& pos_j) {
    /**
     * Calculate the Euclidean distance between two 3D positions.
     *
     * @param pos_i: 3D position of particle i.
     * @param pos_j: 3D position of particle j.
     * @return: Euclidean distance between pos_i and pos_j.
     */
    double dx = pos_i[0] - pos_j[0];
    double dy = pos_i[1] - pos_j[1];
    double dz = pos_i[2] - pos_j[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}