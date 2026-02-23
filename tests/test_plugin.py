import os
import sys
import glob
import subprocess
import pytest
import mdi
import time
import traceback
import multiprocessing


# This class assists in handling any exceptions associated with multiprocessing
class PytestProcess(multiprocessing.Process):
    def __init__(self, *args, **kwargs):
        multiprocessing.Process.__init__(self, *args, **kwargs)
        self._pconn, self._cconn = multiprocessing.Pipe()
        self._exception = None

    def run(self):
        try:
            multiprocessing.Process.run(self)
            self._cconn.send(None)
        except Exception as e:
            tb = traceback.format_exc()
            self._cconn.send((e, tb))

    @property
    def exception(self):
        if self._pconn.poll():
            self._exception = self._pconn.recv()
        return self._exception


# This function runs a test in another process.
# This is important to ensures that no state is shared between different test runs.
def run_test_process(task, args):
    # Launch the process
    p = PytestProcess(target=task, args=args)
    p.start()

    # Wait for the process to complete
    p.join()

    # Handle any errors
    if p.exception:
        error, traceback = p.exception
        print(traceback)
        raise error



###############################################################
# Test <NAME
###############################################################

def test_name(plugindir):
    # This function contains the actual test, which will be run in another process
    def task(plugindir):
        mdi.MDI_Init( "-name pytest -role DRIVER -method LINK -plugin_path " + str(plugindir) )

        # Launch an instance of the engine library
        mpi_world = None
        mdi_comm = mdi.MDI_Open_plugin("ERI",
                        "-mdi \"-name ERI -role ENGINE -method LINK\"",
                        mpi_world)

        # Determine the name of the code
        mdi.MDI_Send_Command("<NAME", mdi_comm)
        name = mdi.MDI_Recv(mdi.MDI_NAME_LENGTH, mdi.MDI_CHAR, mdi_comm)
        assert name == "ERI"

        mdi.MDI_Close_plugin(mdi_comm)

    # Run the task function in another process
    run_test_process(task, (plugindir,))




###############################################################
# Test >NATOMS
###############################################################

def test_natoms(plugindir):
    # This function contains the actual test, which will be run in another process
    def task(plugindir):
        mdi.MDI_Init( "-name pytest -role DRIVER -method LINK -plugin_path " + str(plugindir) )

        # Launch an instance of the engine library
        mpi_world = None
        mdi_comm = mdi.MDI_Open_plugin("ERI",
                        "-mdi \"-name ERI -role ENGINE -method LINK\"",
                        mpi_world)

        # Get the number of atoms
        natoms = 10
        mdi.MDI_Send_Command(">NATOMS", mdi_comm)
        mdi.MDI_Send(natoms, 1, mdi.MDI_INT, mdi_comm)

        mdi.MDI_Close_plugin(mdi_comm)

    # Run the task function in another process
    run_test_process(task, (plugindir,))
