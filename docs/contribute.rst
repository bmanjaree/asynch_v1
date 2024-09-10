.. _contributing-to-asynch:

Contributing to Asynch
======================

You are here to help on Asynch? Awesome, feel welcome and read the following sections in order to know what and how to work on something. If you get stuck at any point you can create a `ticket on GitHub`_.

.. _ticket on GitHub: https://github.com/Iowa-Flood-Center/asynch/issues

Contributing to development
---------------------------

If you want to deep dive and help out with development on Read the Docs, then first get the project installed locally according to the :ref:`Installation` instruction. After that is done we suggest you have a look at tickets in our issue tracker that are labelled `Good First Bug`_. These are meant to be a great way to get a smooth start and won't put you in front of the most complex parts of the system.

If you are up to more challenging tasks with a bigger scope, then there are a set of tickets with a `Enhancement`_ tag. These tickets have a general overview and description of the work required to finish. If you want to start somewhere, this would be a good place to start. That said, these aren't necessarily the easiest tickets. They are simply things that are explained. If you still didn't find something to work on, search for the `Sprintable`_ label. Those tickets are meant to be standalone and can be worked on ad-hoc.

When contributing code, then please follow the standard Contribution Guidelines set forth at `contribution-guide.org`_.

.. _Enhancement: https://github.com/Iowa-Flood-Center/asynch/issues?direction=desc&labels=enhancement&page=1&sort=updated&state=open
.. _Good First Bug: https://github.com/rtfd/readthedocs.org/issues?q=is%3Aopen+is%3Aissue+label%3A%22Good+First+Bug%22
.. _Sprintable: https://github.com/rtfd/readthedocs.org/issues?q=is%3Aopen+is%3Aissue+label%3ASprintable
.. _contribution-guide.org: http://www.contribution-guide.org/#submitting-bugs

Keeping the documentation updated
---------------------------------

Whenever an update in the code adds, changes or removes elements that affect the user experience (i.e.: changes in input forcing formats, existing hlm models, global file format, etc ), it is expected from the developer to perform the respective updates in the *Read the Docs* documentation.

Here are some tips and explanations regarding the documentation process.

Documentation structure
~~~~~~~~~~~~~~~~~~~~~~~

The documentation publicly available at ``http://asynch.readthedocs.io`` is hosted by `Read the Docs <https://readthedocs.org/>`__ in a project that is maintained by the same maintainers of Asynch.

This documentation is written in ``reStructuredText`` formatting language (*.rst files) and uses Sphinx with Read the Docs libraries for compilation purposes. All the files are relevant for the documentation are located in the ``doc/`` folder.

Working locally
~~~~~~~~~~~~~~~

A desktop environment can be set up to compile the documentation locally. This approach is good for debugging the result before publishing it publicly.

For setting up the environment, you will have to install Read the Docs on your machine in a separated virtual Python environment.

Ensure **git 1.5** or higher, **Python 3.6** or higher, and both Python **virtualenv** and **virtualenvwrapper** are installed. If using Windows, ensure Python **virtualenvwrapper-win** and `cmder <http://cmder.net/>`__ are also installed (use **cmder** to perform the following command operations).

Create a virtual environment for compiling Read the Docs documentations:
  
.. code-block:: sh

   mkvirtualenv readthedocs

If your command session does not activate the newly created virtual environment automatically, activate it:

.. code-block:: sh

   workon readthedocs
   
Navigate to this virtualenv directory:

.. code-block:: sh

   cd [USER_HOME]\Env\readthedocs\

Checkout Read the Docs into a new folder and enter there:

.. code-block:: sh

   mkdir checkouts
   cd checkouts
   git clone https://github.com/rtfd/readthedocs.org.git
   cd readthedocs.org\

Install all requirements for Read The Docs:
  
.. code-block:: sh

   pip install -r requirements.txt
  
or:

.. code-block:: sh

   python -m pip install -r requirements.txt
   
Navigate to the ``docs`` directory of the local clone of Asynch repository:

.. code-block:: sh

   cd [ASYNCH]\docs\
   
Perform the changes you want in the ``.srt`` files within this folder. After that, compile using the command:

.. code-block:: sh

   make html
   
.. note::

   1-) When compiling, ensure you are still working on the readthedocs virtual Python environment;
   
   2-) When compiling, some Python packages may be required. Be prepared to perform pip installs;
   
   3-) The command 'make' also works on Windows when runned within cmder.
   
Access the results opening the file ``[ASYNCH]/docs/.build/html/index.html`` with a web browser.

This mini tutorial was adapted from `here <http://docs.readthedocs.io/en/latest/install.html>`__.
     
Publishing
~~~~~~~~~~

Everytime a git ``push`` or ``pull request`` is performed into the ``master`` or ``develop`` branches in the Git Hub repository, or on a branch that creates a new ``tag``, the Read the Docs server reads, compiles and publishes the documentation online.

This connection between Git Hub and Read the Docs is established through the so called *webhooks*. The official Asynch Git Hub account has a webhook that triggers the compiling steps from the Read the Docs server.

Managing releases
-----------------

Once you are happy with your changes in the ``develop`` branch and ran a couple of test simulations, here is the procedure to release a new version ``x.y.z`` (e.g. ``1.5.0``):

Branch
~~~~~~

Create a branch for the release following the ``release-x.y.z`` naming scheme and `semantic versionning <http://semver.org/>`__ rules :

.. code-block:: sh

  git branch release-x.y.z

Edit
~~~~

Edit the release notes (``doc/release_notes.rst``).

Edit ``configure.ac`` to bump the version number:

.. code::

  AC_INIT([asynch], [x.y.z], [samuel-debionne@uiowa.edu])

Commit your changes.

.. code-block:: sh

  git add configure.ac doc/release_notes.rst
  git commit -m "Bump version number to x.y.z"
  git push

Generate the tarball
~~~~~~~~~~~~~~~~~~~~

In a new empty folder, run the following commands to clone the repository, generate the configure script and the tarball.

.. code-block:: sh

  git clone https://github.com/Iowa-Flood-Center/asynch.git
  git checkout release-x.y.z
  autoreconf -i
  mkdir build && cd build
  export TAR_OPTIONS="--owner=0 --group=0 --numeric-owner"
  ../configure
  make dist

That should generate a ``release-x.y.z.tar.gz`` that needs to be tested.

Test the tarball
~~~~~~~~~~~~~~~~

In a new empty folder, follow  the instructions in :ref:`Installing the package`:

.. code-block:: sh

  tar xf release-x.y.z.tar.gz
  cd release-x.y.z
  mkdir build && cd build
  ../configure CFLAGS="-O2 -DNDEBUG"
  make
  make check
  make install

Adjust the release branch if there is any problem with the build (e.g. missing header file).

Release on Github
~~~~~~~~~~~~~~~~~

Merge the release branch ``release-x.y.z`` to ``master``. The easiest way is to submit a new Pull Request. The *base* branch should be ``Iowa-Flood-Center/asynch`` / ``master`` and the *compare* branch ``Iowa-Flood-Center/asynch`` / ``release-x.y.z``.


Review your Pull Request, or better let someone else do the review. If everything looks good, and you have `Travis CI <https://travis-ci.org/Iowa-Flood-Center/asynch>`_'s blessing, do a *"Merge and Squash"*.

You can safely delete the release branch at this point.

Click on *"Draft a new release"* in `Releases <https://github.com/Iowa-Flood-Center/asynch/releases>`_:

=============== ===============
Field           Value
=============== ===============
Tag version     vx.y.z (v1.5.0)
Release title   Pick a city in Iowa
Description     A short version of the release notes
=============== ===============

Attach the tarball that was generated in the previous step. This is usefull because the tarball does not require the target computer to have autotools installed.

Ready? *"Publish Release"*! Every followers of the repo get notified of the new version. Good job!
