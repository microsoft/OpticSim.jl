# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using OpticSim # this will only work after the main build steps are completed
@info "Running representative workload"
# add stuff here #
Examples.autodrawrays()
Examples.hexapolarspotdiagramexample()
@info "Finished running representative workload"
