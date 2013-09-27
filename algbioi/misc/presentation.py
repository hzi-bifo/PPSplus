import multiprocessing

def _test():

    print("Number of processors: %s" % multiprocessing.cpu_count())

    import sys
    print("Platform: %s" % sys.platform)

    import platform

    print platform.machine()

    print platform.version()

    print platform.platform()

    print platform.uname()

    print platform.system()

    print platform.processor()

    import os
    print os.uname()

    print 'uname:', platform.uname()

    print
    print 'system   :', platform.system()
    print 'node     :', platform.node()
    print 'release  :', platform.release()
    print 'version  :', platform.version()
    print 'machine  :', platform.machine()
    print 'processor:', platform.processor()

    print '------------'
    print platform.architecture()
    print platform.python_build()
    print platform.python_compiler()
    print platform.python_version()
    print platform.python_implementation()
    print platform.uname()
    print platform.linux_distribution()

    print '------------'
    print os.getcwd()



if __name__ == "__main__":
    _test()