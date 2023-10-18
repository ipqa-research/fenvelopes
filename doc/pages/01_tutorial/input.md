---
title: Setting up the config file
---

The input files are based on the `namelist` format provided intrisicly in the
Fortran standard. An usual input file is structured like:

```fortran
! input.nml
!
! Namelist based input file
! =========================
!
! Units:
!  - Pressure: bar
!  - Temperature: K
!  - Volume: L
! =========================


&nml_setup
    ! General settings
    nc=5,                ! Number of components
    model="PR78",        ! SRK PR76 PR78
    mixrule="ClassicVdW" ! only ClassicVdW for now
/

&nml_composition
    names="PC1" "PC2" "PC3" "PC4" "H2O"
    spec="critical", ! critical or parameters specification
    z=0.15 0.10 0.10 0.15 0.50
/

&nml_classicvdw ! Classic VdW mixing rule parameters
    ! kij matrix
    kij(1, :)=0      0      0      0      0.7192
    kij(2, :)=0      0      0      0      0.4598
    kij(3, :)=0      0      0      0      0.2673
    kij(4, :)=0      0      0      0      0.2417
    kij(5, :)=0.7192 0.4598 0.2673 0.2417 0
    
    ! lij matrix
    lij(:, :) = 0
/

&nml_critical
    ! Critical constants
    
    ! Critical Temperature
    tc=305.586 638.889 788.889 838.889 647.3
    
    ! Critical Pressure
    pc=48.82 19.65 10.2 7.72 220.89

    ! Acentric Factor
    w=0.098 0.535 0.891 1.085 0.344
/

&nml_px 
    ! Px envelopes relevant info
    ! Temperature
    T=350.0
    
    ! Initial composition, ussualy the same as the main fluid.
    z_0=0.15 0.10 0.10 0.15 0.50 
    
    ! Injection fluid composition
    z_injection=1 0 0 0 0

    ! Which kind of injection to realize
    injection_case="displace" ! [dilute|displace]
/
```