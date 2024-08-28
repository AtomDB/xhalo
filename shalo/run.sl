#!/usr/bin/env slsh
require("shalo");
require("interpol");

define slsh_main ()
{
    variable NH = 1e22;
    variable E1 = 2.0;
    variable dustmodel = 1;
    variable halomodel = 1;
    variable SrcRad = 0.0;
    variable th = [10,20,30,40,50];
    variable Isca = shalo (NH, E1, dustmodel, halomodel, SrcRad, th);
    print(Isca);
}