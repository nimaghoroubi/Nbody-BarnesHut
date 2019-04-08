%%%%%%%%%% Read me if you want to live! %%%%%%%%%%%%%%%%%%%%%%%%
///demangling
    int N            = atoi(argv[1]);
    char* filename   =     (argv[2]);
    int step_num     = atoi(argv[3]);
    double dt        = atof(argv[4]);
    double theta_max = atof(argv[5]);
    graphics     = atoi(argv[6]);
    
    
    
    this means if you want to run the program, you do it like this ./a.out 1000 ellipse_N_01000.gal 200 1e-5 0.25 1
    this means you are using barnes hut method with theta 0.25 (https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation) 
    so, a.out is executable, 1000 is number of bodies, ellipse_N_01000.gal is file name (i included for your fun!) 200 is the time steps
    1e-5 is time step, 0.25 is theta in barnes hut, 1 means graphics are on! (for more fun even, yay!)
