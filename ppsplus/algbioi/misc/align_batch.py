#!/usr/bin/env python

"""
    Copyright (C) 2014  Ivan Gregor

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Note that we could have written some parts of this code in a nicer way,
    but didn't have time. Be careful when reusing the source code.


    To computed alignments of many FASTA files in parallel.

"""
import os
import sys
import multiprocessing as mp
import subprocess
import tempfile

# SEE BELOW FOR MAIN !!!

# Taken from algbioi.com.parallel.py !!!


class TaskThread():
    def __init__(self, fun, args):
        """
            Defines one function and its arguments to be executed in one thread.

            @param fun: a function to be executed
            @type fun: function
            @param args: arguments of the function
            @type args: tuple
        """
        self.fun = fun
        self.args = args


class TaskCmd():
    def __init__(self, cmd, cwd='.', stdin=None, stdout=None, stderr=None):
        """
            Defines one task to be executed as a command line command.

            @param cmd: command to be executed on a command line
            @param cwd: current working directory in which the task will be executed
            @param stdin: process standard input
            @param stdout: process standard output
            @param stderr: process standard err
        """
        self.cmd = cmd
        self.cwd = cwd
        self.stdin = stdin
        self.stdout = stdout
        self.stderr = stderr


def runThreadParallel(threadTaskList, maxThreads=mp.cpu_count()):
    """
        Execute several functions (threads, processes) in parallel.

        @type threadTaskList: list of TaskThread
        @param maxThreads: maximum number of tasks that will be run in parallel at the same time
        @return: a list of respective return values
    """
    assert isinstance(threadTaskList, list)
    assert isinstance(maxThreads, int)

    # creates a pool of workers, add all tasks to the pool
    pool = mp.Pool(processes=maxThreads)
    taskHandlerList = []
    for task in threadTaskList:
        assert isinstance(task, TaskThread)
        taskHandlerList.append(pool.apply_async(task.fun, task.args))

    # finish all tasks
    pool.close()
    pool.join()

    # retrieve the return values
    retValList = []
    for taskHandler in taskHandlerList:
        taskHandler.wait()
        assert taskHandler.successful()
        retValList.append(taskHandler.get())

    return retValList


def _runCmd(taskCmd, stdInErrLock=None):
    """
        Executes a command line task.

        @type taskCmd: TaskCmd
        @param stdInErrLock: acquiring the lock enables writing to the stdout and stderr (if not None)
        @type stdInErrLock: multiprocessing.Lock

        @return: a tuple (process, TaskCmd)
    """
    # setting up stdin and stdout (to buffer the output)
    if taskCmd.stdout is None and stdInErrLock is not None:
        stdout = tempfile.TemporaryFile(mode='w+r')
        stdoutP = stdout
    else:
        stdout = None
        stdoutP = taskCmd.stdout

    if taskCmd.stderr is None and stdInErrLock is not None:
        stderr = tempfile.TemporaryFile(mode='w+r')
        stderrP = stderr
    else:
        stderr = None
        stderrP = taskCmd.stderr

    # running the command line task
    try:
        process = subprocess.Popen(taskCmd.cmd, shell=True, bufsize=-1, cwd=taskCmd.cwd, stdin=taskCmd.stdin,
                                   stdout=stdoutP, stderr=stderrP)
        process.wait()
    finally:
        # exclusive writing to the stdin or stderr (empty the buffers containing stdin or stdout of the run)
        if stdout is not None or stderr is not None:
            stdInErrLock.acquire()
            if stdout is not None:
                stdout.flush()
                stdout.seek(0)
                sys.stdout.write(stdout.read())
                sys.stdout.flush()
                stdout.close()
            if stderr is not None:
                stderr.flush()
                stderr.seek(0)
                sys.stderr.write(stderr.read())
                sys.stderr.flush()
                stderr.close()
            stdInErrLock.release()

    return (process, taskCmd)


def runCmdParallel(cmdTaskList, maxProc=mp.cpu_count(), stdInErrLock=mp.Manager().Lock()):
    """
        Run several command line commands in parallel.

        @attention: use the Manager to get the lock as in this function definition !!!

        @param cmdTaskList: list of command line tasks
        @type cmdTaskList: list of TaskCmd
        @param maxProc: maximum number of tasks that will be run in parallel at the same time
        @param stdInErrLock: acquiring the lock enables writing to the stdout and stderr

        @return: list of failed commands, dictionary (cmd, task process)
    """
    assert isinstance(cmdTaskList, list)
    assert isinstance(maxProc, int)

    threadTaskList = []
    for cmdTask in cmdTaskList:
        assert isinstance(cmdTask, TaskCmd)

        threadTaskList.append(TaskThread(_runCmd, (cmdTask, stdInErrLock)))

    returnValueList = runThreadParallel(threadTaskList, maxProc)

    failList = []
    for process, task in returnValueList:
        if process.returncode != 0:
            failList.append(dict(process=process, task=task))
    if len(failList) > 0:
        return failList
    else:
        return None


def reportFailedCmd(failList):
    """
        Report on failed commands.
    """
    if failList is not None:
        assert isinstance(failList, list)
        msgList = []
        for task in failList:
            assert isinstance(task, dict)
            msg = 'Task failed with return code: %s, task: %s' % (task['process'].returncode, task['task'].cmd)
            msgList.append(msg)
            sys.stderr.write(msg)
        sys.stderr.flush()
        return msgList
    else:
        return None


# MAIN !!!


def getAlignments(inDir, outDir, muscleBinary, maxCPU):
    """
        Build a multiple sequence alignments for each file in the input directory.

    """
    assert os.path.isfile(muscleBinary), 'Binnary file does not exist: %s' % muscleBinary
    taskList = []
    for fileName in os.listdir(inDir):
        cmd = '%s -in %s -out %s' % (muscleBinary, os.path.join(inDir, fileName), os.path.join(outDir, fileName))
        taskList.append(TaskCmd(cmd, cwd=outDir))

    reportFailedCmd(runCmdParallel(taskList, maxProc=maxCPU))


if __name__ == "__main__":
    # muscleBinary = '/net/metagenomics/projects/PPSmg/hsim/muscle3.8.31_i86linux64'
    muscleBinary = '/Users/ivan/Documents/work/tools/muscle/muscle3.8.31_i86darwin64'
    # inDir = '/net/metagenomics/projects/PPSmg/hsim/phylo_acc_all_genes'
    inDir = '/Users/ivan/Documents/nobackup/hsim01/562/a'
    # outDir = '/net/metagenomics/projects/PPSmg/hsim/phylo_acc_all_genes_align_a'
    outDir = '/Users/ivan/Documents/nobackup/hsim01/562/b'
    getAlignments(inDir, outDir, muscleBinary, maxCPU=mp.cpu_count())