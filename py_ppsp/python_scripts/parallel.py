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

    Provides convenient functions to run functions and command line commands in parallel.
"""

import os
import sys
import time
import multiprocessing as mp
import subprocess
import tempfile


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


def runThreadParallel(threadTaskList, maxThreads=mp.cpu_count(), keepRetValues=True):
    """
        Execute several functions (threads, processes) in parallel.

        @type threadTaskList: list[TaskThread]
        @param maxThreads: maximum number of tasks that will be run in parallel at the same time
        @param keepRetValues: keep the return values and return them at the end (else return None)
        @type keepRetValues: bool
        @return: a list of respective return values (or None)
        @rtype: list | None
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
    if keepRetValues:
        retValList = []
    else:
        retValList = None

    for taskHandler in taskHandlerList:
        taskHandler.wait()
        assert taskHandler.successful()
        if keepRetValues:
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
    closeStdout = False
    if taskCmd.stdout is not None and type(taskCmd.stdout) is str:
        taskCmd.stdout = open(taskCmd.stdout, 'w')

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

    if closeStdout:
        taskCmd.stdout.close()

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


def runCmdSerial(cmdTaskList, verbose=False, stopWhenError=True, stdInErrLock=None):
    """
        Run several command line commands one by one.

        @attention: Use the Manager to get the lock (mp.Manager().Lock()) if the lock shared among multiple processes!

        @param cmdTaskList: list of command line tasks
        @type cmdTaskList: list of TaskCmd
        @param stdInErrLock: acquiring the lock enables writing to the stdout and stderr
        @type stdInErrLock: multiprocessing.Lock()
    """
    assert isinstance(cmdTaskList, list)

    counter = 0
    failList = []
    for task in cmdTaskList:
        counter += 1
        if verbose:
            msg = 'Starting "#%s" cmd: %s\n' % (counter, task.cmd)
            if stdInErrLock is not None:
                stdInErrLock.acquire()
            sys.stdout.write(msg)
            sys.stdout.flush()
            if stdInErrLock is not None:
                stdInErrLock.release()

        # run command
        process, taskCmd = _runCmd(task, stdInErrLock)

        if process.returncode != 0:
            failList.append(dict(process=process, task=task))
            if stopWhenError:
                break
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


# Deprecated implementation!


class _ProcHandlerCmd():
    """
        Stores context of one process.

        @deprecated: there is a newer implementation !!!
    """
    def __init__(self, task, timeout=None):
        """
            @type task: TaskCmd
        """
        self.process = subprocess.Popen('exec ' + task.cmd, shell=True, bufsize=-1, cwd=task.cwd)
        self.cmd = task.cmd
        self.runtime = 0.
        self.timeout = timeout

    def incRuntime(self, timeStep):
        self.runtime += timeStep

    def isTimeOut(self):
        if self.timeout is not None and self.runtime > self.timeout:
            return True
        return False

    def getPid(self):
        return self.process.pid


def runCmdParallel0(cmdTaskList, maxProc=mp.cpu_count(), timeout=None, timeStep=0.3):
    """
        Run several command line commands in parallel.

        @deprecated: there is a newer implementation !!!
        @see: runCmdParallel
        @warning: use just for testing purposes when making use of the timeout option

        @param cmdTaskList: list of command line tasks
        @type cmdTaskList: list of TaskCmd
        @param maxProc: maximum number of tasks that will be run in parallel at the same time
        @param timeout: after this number of seconds, the process will be killed, (None if no timeout set)
        @param timeStep: time interval in which processes will be actively checked whether they are running
        @return: list of failed commands, tuple (command, return code, runtime, pid)
    """
    counter = 0
    failList = []
    cmdInGen = True
    cmdGen = iter(cmdTaskList)  # generator of commands
    procArray = {}
    for i in range(maxProc):
        procArray[i] = None

    # loop until all processes finish (or are killed)
    while True:

        # run commands
        if cmdInGen:
            for i in range(maxProc):
                if procArray[i] is None:
                    try:
                        task = cmdGen.next()
                        procArray[i] = _ProcHandlerCmd(task, timeout)  # run process
                        counter += 1
                        print('Running "%s" cmd: %s' % (counter, task.cmd))
                    except StopIteration:
                        cmdInGen = False  # there are no processes to be run

        # sleep for a while
        time.sleep(timeStep)

        # check for finished processes and processes passed timeout
        for i in range(maxProc):
            ph = procArray[i]
            if ph is not None:  # there is a process in the slot
                if ph.process.poll() is None:
                    # process running
                    ph.incRuntime(timeStep)
                    if ph.isTimeOut():
                        # process over time, kill it!
                        ph.process.kill()
                        print("Process (%s): %s killed! (after %ss)" % (ph.getPid(), ph.cmd, ph.runtime))
                        failList.append((ph.cmd, 9, ph.runtime, ph.getPid()))
                        procArray[i] = None  # free slot
                else:
                    # process finished
                    ph.process.wait()
                    if ph.process.returncode != 0:
                        print('Process(%s): "%s" ended with return code: "%s' % (ph.getPid(), ph.cmd,
                                                                                 ph.process.returncode))
                        failList.append((ph.cmd, ph.process.returncode, ph.runtime, ph.getPid()))
                    procArray[i] = None  # free slot

        # finish if no process is running and there is not process to be run
        if len(set(procArray.values())) == 1 and not cmdInGen:
            break

    return failList


# TESTING !!!

def _testCmd(parallel=True):
    print('Start: Test: runCmdParallel')
    inDir = '/Users/ivan/Documents/nobackup/hsim01/562/a'
    outDir = '/Users/ivan/Documents/nobackup/hsim01/562/b'
    MUSCLE_BINARY = '/Users/ivan/Documents/work/tools/muscle/muscle3.8.31_i86darwin64'
    assert os.path.isfile(MUSCLE_BINARY), 'Binnary file does not exist: %s' % MUSCLE_BINARY
    cmdListA = []
    for fileName in os.listdir(inDir):
        cmd = '%s -in %s -out %s' % (MUSCLE_BINARY, os.path.join(inDir, fileName), os.path.join(outDir, fileName))
        # print cmd
        cmdListA.append(TaskCmd(cmd, outDir))
        # break

    if parallel:
        failList = runCmdParallel(cmdListA)
    else:
        lock = mp.Lock()
        failList = runCmdSerial(cmdListA, stdInErrLock=lock)
    reportFailedCmd(failList)
    print('Stop: Test: runCmdParallel')


def _f(a, t, n, c):
    for i in range(n):
        print(a)
        cc = 333
        for j in range(int(10000000*t)):
            c /= cc
    return t, c


def _testThread():
    print('Start: Test: RunThreadParallel')
    r = runThreadParallel([TaskThread(_f, ('thread 1', 0.5, 5, 48749394857384234987)),
                           TaskThread(_f, ('thread 2', 0.7, 6, 57395769304867332349)),
                           TaskThread(_f, ('thread 3', 0.8, 7, 87263485768798234987)),
                           TaskThread(_f, ('thread 4', 0.9, 8, 38573947573957684485)),
                           TaskThread(_f, ('thread 5', 0.9, 8, 38573947573957684485)),
                           TaskThread(_f, ('thread 6', 1.0, 8, 38573947573957684485)),
                           TaskThread(_f, ('thread 7', 1.1, 8, 38573947573957684485))])
    print(r)
    print('Stop: Test: RunThreadParallel')


def _testMisc():
    cmd = 'echo "a; echo "b" >&2'
    lock = mp.Manager().Lock()

    # print runCmdSerial([TaskCmd(cmd)], stdInErrLock=lock, verbose=True)
    t = TaskThread(_runCmd, (TaskCmd(cmd), lock))
    print runThreadParallel([t])


# if __name__ == "__main__":
#     pass
    # reportFailedCmd(runCmdSerial([TaskCmd('pythonn -version')], verbose=True))
    # _testThread()
    # _testCmd()
    # _testMisc()
