# include "encabezados.h"

# define GRO_FLAG 1
# define HOST 0
# define FACTOR 2.0

//===================================== 2025-03-11 ===================================

bool between(const float3 rr1, const float3 rr2, const float3 *rr_sml_vec, parametros pars){

    //unitary vector
    float3 u, rr, dr1, dr2, uxrr;
    float u_abs, dr1Norm, dr2Norm, dr_line;
    float max_distance_big = 0.5 * pars.sigma_big + FACTOR * pars.sigma_sml;
    float max_distance_line = 0.35 * pars.sigma_big;
    float side, side_inv;
    side = pars.side;
    side_inv = 1 / side;

    u.x = rr1.x -rr2.x;
    u.y = rr1.y -rr2.y;
    u.z = rr1.z -rr2.z;
    // periodic conditions
    // periodic boundary conditions
    u.x -= side*floor(side_inv*u.x + 0.5f);
    u.y -= side*floor(side_inv*u.y + 0.5f);
    u.z -= side*floor(side_inv*u.z + 0.5f);

    u_abs = sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
    u.x /= u_abs;
    u.y /= u_abs;
    u.z /= u_abs;

    for (int i = 0; i < pars.ngrain_sml; ++i){
        rr = rr_sml_vec[i];
        dr1.x = rr.x - rr1.x;
        dr1.y = rr.y - rr1.y;
        dr1.z = rr.z - rr1.z;
        // periodic boundary conditions
        dr1.x -= side*floor(side_inv*dr1.x + 0.5f);
        dr1.y -= side*floor(side_inv*dr1.y + 0.5f);
        dr1.z -= side*floor(side_inv*dr1.z + 0.5f);
        dr1Norm = sqrt(dr1.x * dr1.x + dr1.y * dr1.y + dr1.z * dr1.z);
        if (dr1Norm > max_distance_big) continue;
        dr2.x = rr.x - rr2.x;
        dr2.y = rr.y - rr2.y;
        dr2.z = rr.z - rr2.z;
        // periodic boundary conditions
        dr2.x -= side*floor(side_inv*dr2.x + 0.5f);
        dr2.y -= side*floor(side_inv*dr2.y + 0.5f);
        dr2.z -= side*floor(side_inv*dr2.z + 0.5f);
        dr2Norm = sqrt(dr2.x * dr2.x + dr2.y * dr2.y + dr2.z * dr2.z);
        if (dr2Norm > max_distance_big) continue;
        // calculate distance to the line
        uxrr.x = dr1.y * u.z - dr1.z * u.y;
        uxrr.y = -dr1.x * u.z + dr1.z * u.x;
        uxrr.z = dr1.x * u.y - dr1.y * u.x;
        dr_line = sqrt(uxrr.x * uxrr.x + uxrr.y * uxrr.y + uxrr.z * uxrr.z);
        if (dr_line > max_distance_line) continue;
        return true;
    }
    return false;
}

void color_selection(char *color, const bool some_near, const int color_id, int &blue_count, int &green_count, int &red_count){
    if (some_near){
        if (color_id > 0){
            sprintf(color, "Red");
            red_count++;
        }
        else{
            sprintf(color, "Green");
            green_count++;
        }
    }
    else{
        sprintf(color,"Red");
        blue_count++;
    }
}

int main() {
    float3 *rr_big_vec, *rr_big_raw_vec, *rr_big_ini_vec, *vv_big_vec, *ff_big_vec,
            *rr_sml_vec, *rr_sml_raw_vec, *rr_sml_ini_vec, *vv_sml_vec, *ff_sml_vec;
    float3 rr, rrn, rrm, vv, drr, zero;
    float *vir_big_vec, *vir_sml_vec, *pot_big_vec, *pot_sml_vec,
            energs[2];
    float *gder_bb_vec, *gder_bs_vec, *gder_sb_vec, *gder_ss_vec;
    int *nocup_big_vec, *cell_big_vec, *nocup_sml_vec, *cell_sml_vec;
    double run_time, trans_time, dt;
    float side, side_inv, cutoff, bin_size_gder,
            gder_bb, gder_bs, gder_sb, gder_ss, shell_vol, vol_free, msd_big,
            msd_sml, dist, aux, sigma, sigma_big, sigma_sml, mass_big, mass_sml,
            range_gder, xngrain_big, xngrain_sml, xngrain_tot, big_z,
            big_z_ave, virial, ene_pot_big, ene_pot_sml, xnb, T0, Tf, t0, tf, time,
            cell_side_big, phi_big, phi_sml, cell_side_sml, volume, diameter;
    float ***gders, *energy_kb, *energy_ks, *energy_ub, *energy_us, *energy_temp, *msd_b, *msd_s,
            *time_energy, *time_msd;
    int ngrain_tot, ngrain_big, ngrain_sml, ni, niter, ntrans, ngap,
            nsamples, mm, nn, ncell_big, ncell_sml, ncell_big3, ncell_sml3, ntot_big,
            ntot_sml, ii, jj, nbins_gder, nb, ntags_big,
            ntags_sml, nocup, nocup_big_max, nocup_sml_max, counter, n_rescalings_per_gap, ngap_rescaling, n_configs;
    int NH, NB_BIG, NB_SML;
# if !HOST
    int NB_CELL_BIG3, NB_CELL_SML3, NB_NTOT_BIG, NB_NTOT_SML;
# endif
    parametros pars;
    char renglon[200], infile[80], gder_fn[80], energy_fn[80], msd_fn[80], snapshots_fn[80];
    FILE *fp_bitac, *fp_snaps, *fp_energ, *fp_force, *fp_gder, *fp_msd, *fp_in, *fp_out,
            *fp_press, *fp_data, *fp_colors;
    fp_data = fopen("wca.data", "r");
    if (fp_data == nullptr) fp_data = stdin;

    printf("temp_settings ?\n");
    fgets(renglon, sizeof(renglon), fp_data);
    sscanf(renglon, "%f  %f  %f  %f  %d", &T0, &Tf, &t0, &tf, &ngap_rescaling);

    printf("dt ?\n");
    fgets(renglon, sizeof(renglon), fp_data);
    sscanf(renglon, "%lf", &dt);

    printf("trans_time, run_time, nsamples (taken in run)\n");
    fgets(renglon, sizeof(renglon), fp_data);
    sscanf(renglon, "%lf %lf %d", &trans_time, &run_time, &nsamples);

    printf("ntags: big, sml ?\n");
    fgets(renglon, sizeof(renglon), fp_data);
    sscanf(renglon, "%d %d", &ntags_big, &ntags_sml);

    printf("first infile ?\n");
    fgets(renglon, sizeof(renglon), fp_data);
    sscanf(renglon, "%s", infile);

    printf("n_configs ?\n");
    fgets(renglon, sizeof(renglon), fp_data);
    sscanf(renglon, "%d", &n_configs);

    printf("nbins_gder, range_gder ?\n");
    fgets(renglon, sizeof(renglon), fp_data);
    sscanf(renglon, "%d %f", &nbins_gder, &range_gder);

    fclose(fp_data);
    printf("Starts simulation_3D_2SP Device\n");

    // calcula y almacena constantes

    pars.dt = dt;
    pars.temp_set = T0;
    pars.nbins_gder = nbins_gder;
    pars.range_gder = range_gder;
    pars.ntags_big = ntags_big;
    pars.ntags_sml = ntags_sml;

    bin_size_gder = range_gder / ((float) nbins_gder);
    pars.bin_size_gder = bin_size_gder;

    zero.x = zero.y = zero.z = 0.0;

    //Initialize memory to accumulate results

    //Values to calculate size of arrays
    niter = (int) (run_time / dt);
    ngap = niter / nsamples;
    n_rescalings_per_gap = ngap / ngap_rescaling;

    //g_ij(r)
    gders = (float ***) malloc(nsamples * sizeof(float **));
    for (int i = 0; i < nsamples; ++i) {
        gders[i] = (float **) malloc(8 * sizeof(float *));
        for (int j = 0; j < 8; ++j) {
            gders[i][j] = (float *) malloc(nbins_gder * sizeof(float));
            for (int k = 0; k < nbins_gder; ++k) gders[i][j][k] = 0.0;
        }
    }

    //energy and msd
    int n_data_energy = nsamples * n_rescalings_per_gap;
    time_energy = (float *) malloc(n_data_energy * sizeof(float));
    energy_kb = (float *) calloc(n_data_energy, sizeof(float));
    energy_ks = (float *) calloc(n_data_energy, sizeof(float));
    energy_ub = (float *) calloc(n_data_energy, sizeof(float));
    energy_us = (float *) calloc(n_data_energy, sizeof(float));
    energy_temp = (float *) calloc(n_data_energy, sizeof(float));

    time_msd = (float *) malloc(nsamples * sizeof(float));
    msd_b = (float *) calloc(nsamples, sizeof(float));
    msd_s = (float *) calloc(nsamples, sizeof(float));

#if GRO_FLAG
    fp_colors = fopen("color_statistics.out", "w");
#endif
    ////////////////////////////Starts loop//////////////////////////////////////////////////////////////////////

    for (int i_config = 1; i_config <= n_configs; ++i_config) {

        counter = 0;
        int energy_counter = 0;

        // lee datos iniciales
        sprintf(infile, "../configs/init_config_%d", i_config);
        fp_in = fopen(infile, "r");
        if (fp_in == nullptr) {
            printf("Verify file path");
            fflush(stdout);
            exit(-2);
        }

        fgets(renglon, sizeof(renglon), fp_in);
        sscanf(renglon, "%f %f", &phi_big, &phi_sml);
        fgets(renglon, sizeof(renglon), fp_in);
        sscanf(renglon, "%f", &side);
        fgets(renglon, sizeof(renglon), fp_in);
        sscanf(renglon, "%d %d %d", &ngrain_big, &ngrain_sml, &ngrain_tot);
        fgets(renglon, sizeof(renglon), fp_in);
        sscanf(renglon, "%f %f", &sigma_big, &sigma_sml);
        fgets(renglon, sizeof(renglon), fp_in);
        sscanf(renglon, "%f %f", &mass_big, &mass_sml);
        fgets(renglon, sizeof(renglon), fp_in);

        // almacena datos iniciales

        pars.side = side;
        pars.ngrain_big = ngrain_big;
        pars.ngrain_sml = ngrain_sml;
        pars.sigma_big = sigma_big;
        pars.sigma_sml = sigma_sml;
        pars.mass_big = mass_big;
        pars.mass_sml = mass_sml;
        xngrain_big = (float) ngrain_big;
        xngrain_sml = (float) ngrain_sml;
        xngrain_tot = (float) ngrain_tot;
        volume = side * side * side;
        side_inv = 1 / side;

        fp_bitac = fopen("bitacora", "w");

        // cell parameters

        cutoff = sigma_big;
        ncell_big = (int) (side / (1.05 * cutoff));
        cell_side_big = side / ((float) ncell_big);
        pars.ncell_big = ncell_big;
        pars.cell_side_big = cell_side_big;

        cutoff = sigma_sml;
        ncell_sml = (int) (side / (1.05 * cutoff));
        cell_side_sml = side / ((float) ncell_sml);
        pars.ncell_sml = ncell_sml;
        pars.cell_side_sml = cell_side_sml;

        fprintf(fp_bitac, "side  %f  ngr b s tot  %d  %d  %d\n", side, ngrain_big,
                ngrain_sml, ngrain_tot);
        fprintf(fp_bitac, "ncell b, s  %d  %d\n", ncell_big, ncell_sml);
        fprintf(fp_bitac, "cell_side b, s  %f  %f\n", cell_side_big, cell_side_sml);
        fprintf(fp_bitac, "ntags b, s  %d  %d\n", ntags_big, ntags_sml);
        fflush(fp_bitac);

        // ranges

        pars.nrange_bb = 1;
        aux = (sigma_big + sigma_sml) / 2.0;
        pars.nrange_bs = 1 + (int) (aux / cell_side_sml);
        pars.nrange_sb = 1;
        pars.nrange_ss = 1;

        fprintf(fp_bitac, "ranges bb, bs  %d  %d\n", pars.nrange_bb, pars.nrange_bs);
        fprintf(fp_bitac, "ranges sb, ss  %d  %d\n", pars.nrange_sb, pars.nrange_ss);

        // parameters

        ncell_big3 = ncell_big * ncell_big * ncell_big;
        ncell_sml3 = ncell_sml * ncell_sml * ncell_sml;
        ntot_big = ncell_big3 * ntags_big;
        ntot_sml = ncell_sml3 * ntags_sml;
        fprintf(fp_bitac, "ncell_big3  %d  ntot_big  %d\n", ncell_big3, ntot_big);
        fprintf(fp_bitac, "ncell_sml3  %d  ntot_sml  %d\n", ncell_sml3, ntot_sml);

        // memory allocation
#if !HOST

        cudaMallocManaged(&rr_big_vec, ngrain_big*sizeof(float3));
        cudaMallocManaged(&rr_big_raw_vec, ngrain_big*sizeof(float3));
        cudaMallocManaged(&rr_big_ini_vec, ngrain_big*sizeof(float3));
        cudaMallocManaged(&vv_big_vec, ngrain_big*sizeof(float3));
        cudaMallocManaged(&ff_big_vec, ngrain_big*sizeof(float3));
        cudaMallocManaged(&vir_big_vec, ngrain_big*sizeof(float));
        cudaMallocManaged(&pot_big_vec, ngrain_big*sizeof(float));

        cudaMallocManaged(&nocup_big_vec, ncell_big3*sizeof(int));
        cudaMallocManaged(&cell_big_vec, ntot_big*sizeof(int));

        cudaMallocManaged(&rr_sml_vec, ngrain_sml*sizeof(float3));
        cudaMallocManaged(&rr_sml_raw_vec, ngrain_sml*sizeof(float3));
        cudaMallocManaged(&rr_sml_ini_vec, ngrain_sml*sizeof(float3));
        cudaMallocManaged(&vv_sml_vec, ngrain_sml*sizeof(float3));
        cudaMallocManaged(&ff_sml_vec, ngrain_sml*sizeof(float3));
        cudaMallocManaged(&vir_sml_vec, ngrain_sml*sizeof(float));
        cudaMallocManaged(&pot_sml_vec, ngrain_sml*sizeof(float));

        cudaMallocManaged(&nocup_sml_vec, ncell_sml3*sizeof(int));
        cudaMallocManaged(&cell_sml_vec, ntot_sml*sizeof(int));
#else

        rr_big_vec = (float3 *) malloc(ngrain_big * sizeof(float3));
        rr_big_raw_vec = (float3 *) malloc(ngrain_big * sizeof(float3));
        rr_big_ini_vec = (float3 *) malloc(ngrain_big * sizeof(float3));
        vv_big_vec = (float3 *) malloc(ngrain_big * sizeof(float3));
        ff_big_vec = (float3 *) malloc(ngrain_big * sizeof(float3));
        vir_big_vec = (float *) malloc(ngrain_big * sizeof(float));
        pot_big_vec = (float *) malloc(ngrain_big * sizeof(float));

        nocup_big_vec = (int *) malloc(ncell_big3 * sizeof(int));
        cell_big_vec = (int *) malloc(ntot_big * sizeof(int));

        rr_sml_vec = (float3 *) malloc(ngrain_sml * sizeof(float3));
        rr_sml_raw_vec = (float3 *) malloc(ngrain_sml * sizeof(float3));
        rr_sml_ini_vec = (float3 *) malloc(ngrain_sml * sizeof(float3));
        vv_sml_vec = (float3 *) malloc(ngrain_sml * sizeof(float3));
        ff_sml_vec = (float3 *) malloc(ngrain_sml * sizeof(float3));
        vir_sml_vec = (float *) malloc(ngrain_sml * sizeof(float));
        pot_sml_vec = (float *) malloc(ngrain_sml * sizeof(float));

        nocup_sml_vec = (int *) malloc(ncell_sml3 * sizeof(int));
        cell_sml_vec = (int *) malloc(ntot_sml * sizeof(int));
#endif

        // blocks and threads

        NH = 256;
        pars.NH = NH;

        NB_BIG = 1 + (ngrain_big - 1) / NH;
        NB_SML = 1 + (ngrain_sml - 1) / NH;
#if !HOST
        NB_CELL_BIG3 = 1 + (ncell_big3 - 1)/NH;
        NB_CELL_SML3 = 1 + (ncell_sml3 - 1)/NH;
        NB_NTOT_BIG = 1 + (ntot_big - 1)/NH;
        NB_NTOT_SML = 1 + (ntot_sml - 1)/NH;
# endif

        // gets memory for gder.
#if !HOST
        cudaMallocManaged(&gder_bb_vec, NB_BIG*nbins_gder*sizeof(float));
        cudaMallocManaged(&gder_ss_vec, NB_SML*nbins_gder*sizeof(float));
        cudaMallocManaged(&gder_bs_vec, NB_BIG*nbins_gder*sizeof(float));
        cudaMallocManaged(&gder_sb_vec, NB_SML*nbins_gder*sizeof(float));
#else
        gder_bb_vec = (float *) malloc(NB_BIG * nbins_gder * sizeof(float));
        gder_ss_vec = (float *) malloc(NB_SML * nbins_gder * sizeof(float));
        gder_bs_vec = (float *) malloc(NB_BIG * nbins_gder * sizeof(float));
        gder_sb_vec = (float *) malloc(NB_SML * nbins_gder * sizeof(float));

#endif
        // positions, velocities

        for (mm = 0; mm < ngrain_tot; mm++) {
            fgets(renglon, sizeof(renglon), fp_in);
            sscanf(renglon, "%d %f %f %f %f %f %f", &nn, &(rr.x), &(rr.y), &(rr.z),
                   &(vv.x), &(vv.y), &(vv.z));

            if (nn != mm) {
                printf("error: mm %d  nn %d no match\n", mm, nn);
                exit(1);
            }
            if (mm < ngrain_big) {
                rr_big_vec[mm] = rr;
                rr_big_raw_vec[mm] = rr;
                vv_big_vec[mm] = vv;
            } else {
                rr_sml_vec[mm - ngrain_big] = rr;
                rr_sml_raw_vec[mm - ngrain_big] = rr;
                vv_sml_vec[mm - ngrain_big] = vv;
            }
        }
        fclose(fp_in);

        // other parameters

        pars.not_finished = 1;
        niter = (int) (run_time / dt);
        ngap = niter / nsamples;
        if (ngap < 1) {
            printf("error: ngap %d\n", ngap);
            exit(1);
        }
        ntrans = (int) (trans_time / dt);
        fprintf(fp_bitac, "niter, ntrans, nsamples, ngap %d %d %d %d\n", niter, ntrans,
                nsamples, ngap);

        //==================================================================================
        // transient and run
        //==================================================================================

        // clear cell vectors and do locate

#if HOST
        set_vec_int_hst(nocup_big_vec, ncell_big3, 0);
        set_vec_int_hst(cell_big_vec, ntot_big, -1);
        cell_locate_hst('b', rr_big_vec, nocup_big_vec, cell_big_vec, pars);
        set_vec_int_hst(nocup_sml_vec, ncell_sml3, 0);
        set_vec_int_hst(cell_sml_vec, ntot_sml, -1);
        cell_locate_hst('s', rr_sml_vec, nocup_sml_vec, cell_sml_vec, pars);
# else
        set_vec_int_dev<<<NB_CELL_BIG3, NH>>>(nocup_big_vec, ncell_big3, 0);
        set_vec_int_dev<<<NB_NTOT_BIG, NH>>>(cell_big_vec, ntot_big, -1);
        cell_locate_dev<<<NB_BIG, NH>>>('b', rr_big_vec, nocup_big_vec, cell_big_vec, pars);
        set_vec_int_dev<<<NB_CELL_SML3, NH>>>(nocup_sml_vec, ncell_sml3, 0);
        set_vec_int_dev<<<NB_NTOT_SML, NH>>>(cell_sml_vec, ntot_sml, -1);
        cell_locate_dev<<<NB_SML, NH>>>('s', rr_sml_vec, nocup_sml_vec, cell_sml_vec, pars);
# endif

        // clear forces

# if HOST
        set_vec_float3_hst(ff_big_vec, ngrain_big, zero);
        set_vec_float_hst(vir_big_vec, ngrain_big, 0.0);
        set_vec_float_hst(pot_big_vec, ngrain_big, 0.0);
        set_vec_float3_hst(ff_sml_vec, ngrain_sml, zero);
        set_vec_float_hst(vir_sml_vec, ngrain_sml, 0.0);
        set_vec_float_hst(pot_sml_vec, ngrain_sml, 0.0);
# else
        set_vec_float3_dev<<<NB_BIG, NH>>>(ff_big_vec, ngrain_big, zero);
        set_vec_float_dev<<<NB_BIG, NH>>>(vir_big_vec, ngrain_big, 0.0);
        set_vec_float_dev<<<NB_BIG, NH>>>(pot_big_vec, ngrain_big, 0.0);
        set_vec_float3_dev<<<NB_SML, NH>>>(ff_sml_vec, ngrain_sml, zero);
        set_vec_float_dev<<<NB_SML, NH>>>(vir_sml_vec, ngrain_sml, 0.0);
        set_vec_float_dev<<<NB_SML, NH>>>(pot_sml_vec, ngrain_sml, 0.0);
# endif

        // forces

# if HOST
        get_forces_same_hst('b', rr_big_vec, ff_big_vec, vir_big_vec, pot_big_vec,
                            nocup_big_vec, cell_big_vec, pars);
        get_forces_same_hst('s', rr_sml_vec, ff_sml_vec, vir_sml_vec, pot_sml_vec,
                            nocup_sml_vec, cell_sml_vec, pars);
        get_forces_diff_hst('b', 's', rr_big_vec, ff_big_vec, rr_sml_vec, vir_big_vec,
                            pot_big_vec, nocup_sml_vec, cell_sml_vec, pars);
        get_forces_diff_hst('s', 'b', rr_sml_vec, ff_sml_vec, rr_big_vec, vir_sml_vec,
                            pot_sml_vec, nocup_big_vec, cell_big_vec, pars);
# else
        get_forces_same_dev<<<NB_BIG, NH>>>('b', rr_big_vec, ff_big_vec, vir_big_vec,
            pot_big_vec, nocup_big_vec, cell_big_vec, pars);
        get_forces_same_dev<<<NB_SML, NH>>>('s', rr_sml_vec, ff_sml_vec, vir_sml_vec,
            pot_sml_vec, nocup_sml_vec, cell_sml_vec, pars);
        get_forces_diff_dev<<<NB_BIG, NH>>>('b', 's', rr_big_vec, ff_big_vec, rr_sml_vec,
            vir_big_vec, pot_big_vec, nocup_sml_vec, cell_sml_vec, pars);
        get_forces_diff_dev<<<NB_SML, NH>>>('s', 'b', rr_sml_vec, ff_sml_vec, rr_big_vec,
            vir_sml_vec, pot_sml_vec, nocup_big_vec, cell_big_vec, pars);
#endif
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////// run ///////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        big_z_ave = 0.0;
        time = -trans_time;
        for (ni = -ntrans; ni < niter; ni++) {
            cudaDeviceSynchronize();

            if (ni == 0) {
                cudaMemcpy(rr_big_raw_vec, rr_big_vec, ngrain_big * sizeof(float3),
                           cudaMemcpyHostToHost);
                cudaMemcpy(rr_big_ini_vec, rr_big_vec, ngrain_big * sizeof(float3),
                           cudaMemcpyHostToHost);
                cudaMemcpy(rr_sml_raw_vec, rr_sml_vec, ngrain_sml * sizeof(float3),
                           cudaMemcpyHostToHost);
                cudaMemcpy(rr_sml_ini_vec, rr_sml_vec, ngrain_sml * sizeof(float3),
                           cudaMemcpyHostToHost);
            }

            // first half velocity verlet. time increment. HAS TO BE HERE

#if HOST
            update_verlet_init_hst('b', rr_big_vec, rr_big_raw_vec, vv_big_vec, ff_big_vec,
                                   pars);
            update_verlet_init_hst('s', rr_sml_vec, rr_sml_raw_vec, vv_sml_vec, ff_sml_vec,
                                   pars);
# else
            update_verlet_init_dev<<<NB_BIG, NH>>>('b', rr_big_vec, rr_big_raw_vec,
                vv_big_vec, ff_big_vec, pars);
            update_verlet_init_dev<<<NB_SML, NH>>>('s', rr_sml_vec, rr_sml_raw_vec,
                vv_sml_vec, ff_sml_vec, pars);

            cudaDeviceSynchronize();
# endif

            //calculate time and temperature
            time = dt * (1.0 + (float) ni);
            pars.temp_set = calculate_temp_linear(T0, Tf, t0, tf, time);

            // clear cell vectors and do cell locate

# if HOST
            set_vec_int_hst(nocup_big_vec, ncell_big3, 0);
            set_vec_int_hst(cell_big_vec, ntot_big, -1);
            cell_locate_hst('b', rr_big_vec, nocup_big_vec, cell_big_vec, pars);
            set_vec_int_hst(nocup_sml_vec, ncell_sml3, 0);
            set_vec_int_hst(cell_sml_vec, ntot_sml, -1);
            cell_locate_hst('s', rr_sml_vec, nocup_sml_vec, cell_sml_vec, pars);
# else
            set_vec_int_dev<<<NB_CELL_BIG3, NH>>>(nocup_big_vec, ncell_big3, 0);
            set_vec_int_dev<<<NB_NTOT_BIG, NH>>>(cell_big_vec, ntot_big, -1);
            cell_locate_dev<<<NB_BIG, NH>>>('b', rr_big_vec, nocup_big_vec,
                cell_big_vec, pars);
            set_vec_int_dev<<<NB_CELL_SML3, NH>>>(nocup_sml_vec, ncell_sml3, 0);
            set_vec_int_dev<<<NB_NTOT_SML, NH>>>(cell_sml_vec, ntot_sml, -1);
            cell_locate_dev<<<NB_SML, NH>>>('s', rr_sml_vec, nocup_sml_vec,
                 cell_sml_vec, pars);
# endif

            // clear forces

# if HOST
            set_vec_float3_hst(ff_big_vec, ngrain_big, zero);
            set_vec_float_hst(vir_big_vec, ngrain_big, 0.0);
            set_vec_float_hst(pot_big_vec, ngrain_big, 0.0);
            set_vec_float3_hst(ff_sml_vec, ngrain_sml, zero);
            set_vec_float_hst(vir_sml_vec, ngrain_sml, 0.0);
            set_vec_float_hst(pot_sml_vec, ngrain_sml, 0.0);
# else
            set_vec_float3_dev<<<NB_BIG, NH>>>(ff_big_vec, ngrain_big, zero);
            set_vec_float_dev<<<NB_BIG, NH>>>(vir_big_vec, ngrain_big, 0.0);
            set_vec_float_dev<<<NB_BIG, NH>>>(pot_big_vec, ngrain_big, 0.0);
            set_vec_float3_dev<<<NB_SML, NH>>>(ff_sml_vec, ngrain_sml, zero);
            set_vec_float_dev<<<NB_SML, NH>>>(vir_sml_vec, ngrain_sml, 0.0);
            set_vec_float_dev<<<NB_SML, NH>>>(pot_sml_vec, ngrain_sml, 0.0);
# endif

            // get new forces

# if HOST
            get_forces_same_hst('b', rr_big_vec, ff_big_vec, vir_big_vec, pot_big_vec,
                                nocup_big_vec, cell_big_vec, pars);
            get_forces_same_hst('s', rr_sml_vec, ff_sml_vec, vir_sml_vec, pot_sml_vec,
                                nocup_sml_vec, cell_sml_vec, pars);
            get_forces_diff_hst('b', 's', rr_big_vec, ff_big_vec, rr_sml_vec, vir_big_vec,
                                pot_big_vec, nocup_sml_vec, cell_sml_vec, pars);
            get_forces_diff_hst('s', 'b', rr_sml_vec, ff_sml_vec, rr_big_vec, vir_sml_vec,
                                pot_sml_vec, nocup_big_vec, cell_big_vec, pars);
# else
            get_forces_same_dev<<<NB_BIG, NH>>>('b', rr_big_vec, ff_big_vec,
                vir_big_vec, pot_big_vec, nocup_big_vec, cell_big_vec, pars);
            get_forces_same_dev<<<NB_SML, NH>>>('s', rr_sml_vec, ff_sml_vec,
                vir_sml_vec, pot_sml_vec, nocup_sml_vec, cell_sml_vec, pars);
            get_forces_diff_dev<<<NB_BIG, NH>>>('b', 's', rr_big_vec, ff_big_vec,
                rr_sml_vec, vir_big_vec, pot_big_vec, nocup_sml_vec, cell_sml_vec, pars);
            get_forces_diff_dev<<<NB_SML, NH>>>('s', 'b', rr_sml_vec, ff_sml_vec,
                rr_big_vec, vir_sml_vec, pot_sml_vec, nocup_big_vec, cell_big_vec, pars);
# endif

            // finish velocity verlet

# if HOST
            update_verlet_finish_hst('b', vv_big_vec, ff_big_vec, pars);
            update_verlet_finish_hst('s', vv_sml_vec, ff_sml_vec, pars);
# else
            update_verlet_finish_dev<<<NB_BIG, NH>>>('b', vv_big_vec, ff_big_vec, pars);
            update_verlet_finish_dev<<<NB_SML, NH>>>('s', vv_sml_vec, ff_sml_vec, pars);
# endif


            cudaDeviceSynchronize();


            // evaluate energies (0.5 because double count)
            if (ni % ngap_rescaling == 0) {

                //impone temperatura via rescaling
                set_temp(vv_big_vec, vv_sml_vec, energs, pars);

                ene_pot_big = ene_pot_sml = 0.0;
                for (mm = 0; mm < ngrain_big; mm++) ene_pot_big += 0.5 * pot_big_vec[mm];
                ene_pot_big /= xngrain_big;
                for (mm = 0; mm < ngrain_sml; mm++) ene_pot_sml += 0.5 * pot_sml_vec[mm];
                ene_pot_sml /= xngrain_sml;

                energy_ks[energy_counter] += energs[1];
                energy_kb[energy_counter] += energs[0];
                energy_us[energy_counter] += ene_pot_sml;
                energy_ub[energy_counter] += ene_pot_big;
                time_energy[energy_counter] = time;
                energy_temp[energy_counter] = pars.temp_set;

                energy_counter++;
            }

            // processing

            if (ni < 0 && ni % ngap == 0) {
                cudaDeviceSynchronize();

                // get nocup_max

                nocup_big_max = nocup_sml_max = 0;
                for (ii = 0; ii < ncell_big3; ii++) {
                    nocup = nocup_big_vec[ii];
                    if (nocup_big_max < nocup) nocup_big_max = nocup;
                }
                for (ii = 0; ii < ncell_sml3; ii++) {
                    nocup = nocup_sml_vec[ii];
                    if (nocup_sml_max < nocup) nocup_sml_max = nocup;
                }
                printf("run ni %d/%d  %.2f%%  --  print %d/%d  --  ocup b s  %d %d\n",
                       ni, niter, 100.0 * ((float) (ni)) / ((float) niter), counter, nsamples,
                       nocup_big_max, nocup_sml_max);
                fflush(stdout);
            }

            // processing

            if (ni >= 0 && ni % ngap == 0) {
                cudaDeviceSynchronize();
                counter++;

                nocup_big_max = nocup_sml_max = 0;
                for (ii = 0; ii < ncell_big3; ii++) {
                    nocup = nocup_big_vec[ii];
                    if (nocup_big_max < nocup) nocup_big_max = nocup;
                }
                for (ii = 0; ii < ncell_sml3; ii++) {
                    nocup = nocup_sml_vec[ii];
                    if (nocup_sml_max < nocup) nocup_sml_max = nocup;
                }
                printf("run ni %d/%d  %.2f%%  --  print %d/%d  --  ocup b s  %d %d -- temp %.3f\n",
                       ni, niter, 100.0 * ((float) (ni)) / ((float) niter), counter, nsamples,
                       nocup_big_max, nocup_sml_max, pars.temp_set);
                fflush(stdout);
                if (nocup_big_max >= ntags_big || nocup_sml_max >= ntags_sml) {
                    printf("ntags reached limit\n");
                    fflush(stdout);
                    exit(1);
                }

                // writing snapshots

#if GRO_FLAG
                bool high_render=true;
                char color[20];
                int blue_count, green_count, red_count;
                sprintf(snapshots_fn, "snapshot_%d.txt",counter);
                fp_snaps = fopen(snapshots_fn, "w");
                fprintf(fp_snaps, "#include \"colors.inc\"\n");
                if (high_render){
                    fprintf(fp_snaps, "global_settings {\n");
                    fprintf(fp_snaps, "\tradiosity{\n");
                    fprintf(fp_snaps, "\t\tpretrace_start 0.08\n");
                    fprintf(fp_snaps, "\t\tpretrace_end   0.01\n");
                    fprintf(fp_snaps, "\t\tcount 150\n");
                    fprintf(fp_snaps, "\t\tnearest_count 10\n");
                    fprintf(fp_snaps, "\t\terror_bound 0.35\n");
                    fprintf(fp_snaps, "\t\trecursion_limit 2\n");
                    fprintf(fp_snaps, "\t\tlow_error_factor 0.5\n");
                    fprintf(fp_snaps, "\t\tgray_threshold 0.0\n");
                    fprintf(fp_snaps, "\t\tminimum_reuse 0.005\n");
                    fprintf(fp_snaps, "\t\tmaximum_reuse 0.2\n");
                    fprintf(fp_snaps, "\t\tbrightness 1\n");
                    fprintf(fp_snaps, "\t\tadc_bailout 0.005\n");
                    fprintf(fp_snaps, "\t}\n}\n");
                }
                fprintf(fp_snaps, "background {White}\n");
                fprintf(fp_snaps, "camera{\n\tangle 90\n\tlocation <9,9,15>\n\tlook_at <9,9,2>\n}\n");
                fprintf(fp_snaps, "light_source{ <25,25,25> color White}\n");

                blue_count = green_count = red_count = 0;
                for (mm = 0; mm < ngrain_big; mm++) {
                    rrm = rr_big_vec[mm];
                    diameter = sigma_big;
                    int color_id = 0;
                    bool some_near = false;
                    for (nn = 0; nn < ngrain_big; nn++) {
                        if (nn == mm) continue;
                        rrn = rr_big_vec[nn];
                        drr.x = rrn.x - rrm.x;
                        drr.y = rrn.y - rrm.y;
                        drr.z = rrn.z - rrm.z;

                        // periodic boundary conditions
                        drr.x -= side*floor(side_inv*drr.x + 0.5f);
                        drr.y -= side*floor(side_inv*drr.y + 0.5f);
                        drr.z -= side*floor(side_inv*drr.z + 0.5f);

                        dist = drr.x * drr.x + drr.y * drr.y + drr.z * drr.z;
                        if (dist > sigma_big + 2.0 * sigma_sml) continue;
                        between(rrm, rrn, rr_sml_vec, pars) ? color_id++ : color_id--;
                        some_near = true;
                    }
                    color_selection(color, some_near, color_id, blue_count, green_count, red_count);
                    if (1.0 < rrm.z && rrm.z < 3.0)
                        fprintf(fp_snaps, "sphere{ <%8.3f,%8.3f,%8.3f>, %8.3f pigment {%s} finish { specular 0.5 }}\n",
                                rrm.x, rrm.y, rrm.z, diameter * .5, color);

                }
                for (mm = 0; mm < ngrain_sml; mm++)
                {
                    rr = rr_sml_vec[mm];
                    diameter = sigma_sml;
                    if ( 1.2<rr.z && rr.z<2.8) fprintf (fp_snaps, "sphere{ <%8.3f%8.3f%8.3f>, %8.3f pigment {Yellow} finish { specular 0.5 }}\n", rr.x, rr.y, rr.z, diameter*.5);

                }

                fprintf(fp_colors, "%f  %d  %d  %d\n",time, red_count, green_count, blue_count);
                fclose(fp_snaps);
# endif

                //calculate gder

                for (nb = 0; nb < NB_BIG * nbins_gder; nb++) gder_bb_vec[nb] = 0.0f;
                for (nb = 0; nb < NB_SML * nbins_gder; nb++) gder_ss_vec[nb] = 0.0f;
                for (nb = 0; nb < NB_BIG * nbins_gder; nb++) gder_bs_vec[nb] = 0.0f;
                for (nb = 0; nb < NB_SML * nbins_gder; nb++) gder_sb_vec[nb] = 0.0f;

# if HOST
                get_gder_hst('b', 'b', rr_big_vec, rr_big_vec, gder_bb_vec, pars);
                get_gder_hst('b', 's', rr_big_vec, rr_sml_vec, gder_bs_vec, pars);
                get_gder_hst('s', 'b', rr_sml_vec, rr_big_vec, gder_sb_vec, pars);
                get_gder_hst('s', 's', rr_sml_vec, rr_sml_vec, gder_ss_vec, pars);
# else
                get_gder_dev<<<NB_BIG, NH>>>('b', 'b', rr_big_vec, rr_big_vec, gder_bb_vec,
                    pars);
                get_gder_dev<<<NB_BIG, NH>>>('b', 's', rr_big_vec, rr_sml_vec, gder_bs_vec,
                    pars);
                get_gder_dev<<<NB_SML, NH>>>('s', 'b', rr_sml_vec, rr_big_vec, gder_sb_vec,
                    pars);
                get_gder_dev<<<NB_SML, NH>>>('s', 's', rr_sml_vec, rr_sml_vec, gder_ss_vec,
                    pars);

                cudaDeviceSynchronize();
# endif

                for (jj = 0; jj < nbins_gder; jj++) {
                    for (nb = 1; nb < NB_BIG; nb++)
                        gder_bb_vec[jj] += gder_bb_vec[jj + nb * nbins_gder];
                    for (nb = 1; nb < NB_SML; nb++)
                        gder_ss_vec[jj] +=
                                gder_ss_vec[jj + nb * nbins_gder];
                    for (nb = 1; nb < NB_SML; nb++)
                        gder_sb_vec[jj] +=
                                gder_sb_vec[jj + nb * nbins_gder];
                    for (nb = 1; nb < NB_BIG; nb++)
                        gder_bs_vec[jj] +=
                                gder_bs_vec[jj + nb * nbins_gder];
                }


                for (nb = 0; nb < nbins_gder; nb++) {
                    xnb = (float) nb;
                    gder_bb = gder_bs = gder_sb = gder_ss = 0.0f;

                    shell_vol = 4.0 * PI * bin_size_gder * bin_size_gder * bin_size_gder *
                                ((1.0 / 3.0) + xnb * xnb + xnb);

                    sigma = sigma_big;
                    vol_free = volume - (4.0 / 3.0) * PI * sigma * sigma * sigma;
                    gder_bb = gder_bb_vec[nb] / (shell_vol);
                    gder_bb /= (xngrain_big * (xngrain_big - 1.0) / vol_free);

                    sigma = 0.5 * (sigma_big + sigma_sml);
                    vol_free = volume - (4.0 / 3.0) * PI * sigma * sigma * sigma;
                    gder_bs = gder_bs_vec[nb] / (shell_vol);
                    gder_sb = gder_sb_vec[nb] / (shell_vol);
                    gder_bs /= (xngrain_big * xngrain_sml / vol_free);
                    gder_sb /= (xngrain_sml * xngrain_big / vol_free);

                    sigma = sigma_sml;
                    vol_free = volume - (4.0 / 3.0) * PI * sigma * sigma * sigma;
                    gder_ss = gder_ss_vec[nb] / (shell_vol);
                    gder_ss /= (xngrain_sml * (xngrain_sml - 1.0) / vol_free);


                    gders[counter - 1][0][nb] += gder_bb;
                    gders[counter - 1][1][nb] += gder_bs;
                    gders[counter - 1][2][nb] += gder_sb;
                    gders[counter - 1][3][nb] += gder_ss;

                    gders[counter - 1][4][nb] += gder_bb * gder_bb;
                    gders[counter - 1][5][nb] += gder_bs * gder_bs;
                    gders[counter - 1][6][nb] += gder_sb * gder_sb;
                    gders[counter - 1][7][nb] += gder_ss * gder_ss;

                }

                // calculate MSD

                msd_big = 0.0;
                for (mm = 0; mm < ngrain_big; mm++) {
                    drr.x = rr_big_raw_vec[mm].x - rr_big_ini_vec[mm].x;
                    drr.y = rr_big_raw_vec[mm].y - rr_big_ini_vec[mm].y;
                    drr.z = rr_big_raw_vec[mm].z - rr_big_ini_vec[mm].z;
                    msd_big += drr.x * drr.x + drr.y * drr.y + drr.z * drr.z;
                }
                msd_big /= xngrain_big;

                msd_sml = 0.0;
                for (mm = 0; mm < ngrain_sml; mm++) {
                    drr.x = rr_sml_raw_vec[mm].x - rr_sml_ini_vec[mm].x;
                    drr.y = rr_sml_raw_vec[mm].y - rr_sml_ini_vec[mm].y;
                    drr.z = rr_sml_raw_vec[mm].z - rr_sml_ini_vec[mm].z;
                    msd_sml += drr.x * drr.x + drr.y * drr.y + drr.z * drr.z;
                }
                msd_sml /= xngrain_sml;


                time_msd[counter - 1] = time;
                msd_b[counter - 1] = msd_big;
                msd_s[counter - 1] = msd_sml;

            }
        }
        cudaDeviceSynchronize();

        // close files, release stuff and finish

        fclose(fp_bitac);

        cudaFree(rr_big_vec);
        cudaFree(rr_big_raw_vec);
        cudaFree(rr_big_ini_vec);
        cudaFree(vv_big_vec);
        cudaFree(ff_big_vec);
        cudaFree(vir_big_vec);
        cudaFree(pot_big_vec);

        cudaFree(nocup_big_vec);
        cudaFree(cell_big_vec);

        cudaFree(rr_sml_vec);
        cudaFree(rr_sml_raw_vec);
        cudaFree(rr_sml_ini_vec);
        cudaFree(vv_sml_vec);
        cudaFree(ff_sml_vec);
        cudaFree(vir_sml_vec);
        cudaFree(pot_sml_vec);

        cudaFree(nocup_sml_vec);
        cudaFree(cell_sml_vec);

        cudaFree(gder_bb_vec);
        cudaFree(gder_bs_vec);
        cudaFree(gder_sb_vec);
        cudaFree(gder_ss_vec);
        cudaDeviceReset();

        printf("Finished step %d of %d\n", i_config, n_configs);
    }
    //promedia y calcula desviaciones estandar
    fflush(stdout);
    for (int i = 0; i < nsamples; ++i) {

        for (int j = 0; j < 8; ++j) for (int k = 0; k < nbins_gder; ++k) gders[i][j][k] /= n_configs;
        for (int k = 0; k < nbins_gder; ++k) {
            gders[i][4][k] = gders[i][4][k] - gders[i][0][k] * gders[i][0][k];   //sigma^2 = <x^2> - <x>^2
            gders[i][4][k] = sqrt(gders[i][4][k]);

            gders[i][5][k] = gders[i][5][k] - gders[i][1][k] * gders[i][1][k];   //sigma^2 = <x^2> - <x>^2
            gders[i][5][k] = sqrt(gders[i][5][k]);

            gders[i][6][k] = gders[i][6][k] - gders[i][2][k] * gders[i][2][k];   //sigma^2 = <x^2> - <x>^2
            gders[i][6][k] = sqrt(gders[i][6][k]);

            gders[i][7][k] = gders[i][7][k] - gders[i][3][k] * gders[i][3][k];   //sigma^2 = <x^2> - <x>^2
            gders[i][7][k] = sqrt(gders[i][7][k]);
        }
    }

    for (int i = 0; i < n_data_energy; ++i) {
        energy_kb[i] /= n_configs;
        energy_ks[i] /= n_configs;
        energy_ub[i] /= n_configs;
        energy_us[i] /= n_configs;
    }


    //imprime resultados en files
    printf("printing files\n");

    for (int i = 0; i < nsamples; ++i) {

        sprintf(gder_fn, "results/gder_%d.out", i);
        fp_gder = fopen(gder_fn, "w");

        for (nb = 0; nb < nbins_gder; nb++) {
            dist = (0.5 + (float) nb) * bin_size_gder;
            fprintf(fp_gder, "%f  %f  %f  %f  %f  %f  %f  %f  %f\n", dist, gders[i][0][nb], gders[i][1][nb],
                    gders[i][2][nb], gders[i][3][nb], gders[i][4][nb], gders[i][5][nb], gders[i][6][nb],
                    gders[i][7][nb]);

        }
        fclose(fp_gder);
    }

    sprintf(energy_fn, "results/energies_averaged.out");
    fp_energ = fopen(energy_fn, "w");

    for (int i = 0; i < n_data_energy; ++i) {
        fprintf(fp_energ, "%f  %f  %f  %f  %f  %f\n", time_energy[i], energy_ks[i], energy_kb[i], energy_ub[i],
                energy_us[i], energy_temp[i]);
    }
    fclose(fp_energ);

    sprintf(msd_fn, "results/msd_averaged.out");
    fp_msd = fopen(msd_fn, "w");

    for (int i = 0; i < nsamples; ++i) {
        fprintf(fp_msd, "%f  %f  %f\n", time_msd[i], msd_b[i], msd_s[i]);
    }
    fclose(fp_msd);

    for (int i = 0; i < nsamples; ++i) {
        for (int j = 0; j < 8; ++j) free(gders[i][j]);
        free(gders[i]);
    }
    free(gders);

    free(energy_ks);
    free(energy_kb);
    free(energy_us);
    free(energy_ub);
    free(msd_b);
    free(msd_s);

    free(energy_temp);
    free(time_energy);
    free(time_msd);

    fclose(fp_colors);

    printf("TERMINADO\n");
    exit(0);
}
