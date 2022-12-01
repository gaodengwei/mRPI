# mRPI

solve mRPI via SOS

======================I. Toolboxes Required:=================================

1.spotless:
URL: https://github.com/spot-toolbox/spotless

Tips: A fork of the Systems Polynomial Optimization Toolbox. Polynomial computation in the FT is supported by this toolbox.

2.mosek:
URL: http://www.mosek.com

Tips: MOSEK is a highly efficient commercial solver of LPs, QPs, SOCPs, SDPs, and MIPs. This is a basic element of the following toolboxes for solving optimization problems. The toolbox "sedumi" also works well for SDPs but is too slower than mosek.

======================II. the four examples:=================================

1. runLinearCT: linear continuous-time system

2. runLinearDT: linear discrete-time system

3. runNonlinearCT: nonlinear continuous-time system

4. runNonlinearDT: nonlinear discrete-time system

5.runSat/runMK: Attitude Stabilization of A Satellite/ state converge via feedback controller

======================III. Main functions=================================

1.mRPI_CT: for the continuous-time systems

2.mRPI_DT: for the discrete-time systems





