.. _development:

********************
*Prost!* Development
********************

This page describes how to work with *Prost!* in a development environment.

=====
Setup
=====

The setup will walk you through retrieving the code, setting up the git flow
tools, and doing a sample commit.

Grab the *Prost!* codebase from GitHub
--------------------------------------

.. code-block:: bash

   # Choose a repo and grab it:
   # public repo via HTTPS
   git clone https://github.com/uoregon-postlethwait/prost.git
   # public repo via SSH
   git clone git@github.com:uoregon-postlethwait/prost.git
   # development repo via HTTPS
   git clone https://github.com/uoregon-postlethwait/prost-dev.git
   # development repo via SSH
   git clone git@github.com:uoregon-postlethwait/prost-dev.git

   # Enter the development directory.
   cd prost

   # Checkout the branch you want to work on.
   git checkout develop

Optional: Configure local *Prost!* repo for public and development repos
------------------------------------------------------------------------

Some *Prost!* development occurs on a private development GitHub repository
(e.g. quick/dirty custom features for collaborators that aren't fit for public
distribution). 

Setup your local git repo to be able to push/pull from both the public GitHub
*Prost!* repo and the private Github *Prost!* development repo like so.  Your
initial ``.git/config`` should look like this:

.. code-block:: ini
   
   [remote "origin"]
        url = https://github.com/uoregon-postlethwait/prost-dev
        fetch = +refs/heads/*:refs/remotes/origin/*
   [branch "master"]
        remote = origin
        merge = refs/heads/master

Change it to look like this:

.. code-block:: ini
   
   [remote "private"]
        url = https://github.com/uoregon-postlethwait/prost-dev
        fetch = +refs/heads/*:refs/remotes/origin/*
   [remote "public"]
        url = https://github.com/uoregon-postlethwait/prost
        fetch = +refs/heads/*:refs/remotes/origin/*
   [branch "master"]
        remote = private
        merge = refs/heads/master

Setup Git Flow (AVH) Tools
--------------------------

We follow the git flow model of development.  For background see:
`Why Aren't You Using Git Flow? <http://jeffkreeftmeijer.com/2010/why-arent-you-using-git-flow/>`_.

We use the AVH forked version of the git flow tools because the original git
flow tools are no longer maintained, and AVH seems like the best replacement
out there at the moment (it's maintained, and adds a few new features).  
See https://github.com/petervanderdoes/gitflow.

We use a mildly forked version of git flow hooks provided by 
https://github.com/jaspernbrouwer/git-flow-hooks, which provides automatic
version bumping for releases and hotfixes (and a few other small but nice
features).  Our forked/patched version
(see https://github.com/jasonsydes/git-flow-hooks) just adds the tiny feature
of also updating ``prost/_version.py`` in addition to the standard updating of the
``VERSION`` file.

To setup the git flow environment:

* Install git flow (the AVH version):
  https://github.com/petervanderdoes/gitflow/wiki.  Use the defaults for all,
  except for *versiontag*, for which you should specify 'v':

  .. code-block:: bash
     
     Version tag prefix? [] v

* Install our (very mildly) patched version of git-flow-hooks:

  .. code-block:: bash

     # cd to your checkout of the prost repository 
     cd /path/to/prost-repo
     
     # cd into the .git directory
     cd .git
     
     # Rename the hooks directory to something else.  The default git hooks
     # directory just has some examples in it.  Safe to rename, or even delete.
     mv hooks hooks.orig
     
     # Clone our patched version of the git flow hooks repo into hooks.
     git clone https://github.com/jasonsydes/git-flow-hooks hooks

     # cd back to the repo:
     cd ..

HotFix
------

Hotfixes are tiny fixes to an existing release.  The hotfix may optionally be
pushed to the public repo as well. Do the following:

.. code-block:: bash

   # Start a hotfix 
   git flow hotfix start
   # Make your changes to the code (i.e. implement the hotfix).
   vim thecode.py
   # Add and commit the code.
   git add thecode.py
   git ci -m 'Made these changes...'
   # Finish the hotfix
   git flow hotfix finish -m
   # Push everything to private repository
   git push private --mirror
   # Optional: If this is to be a public release, push just the master repo and
   #    the tag just created to the public repo.
   git push public master
   git push public v0.7.17
   
Release
-------

TODO: Add instructions on performing a release with git flow.

Feature
-------

TODO: Add instructions on working on large new features with git flow.

========================================================================
Running *Prost!* from within the development directory (not recommended)
========================================================================

You can run prost directly from within the development directory for quick and
dirty work.  However, we don't recommend this for several reasons:

- Mixing code and data gets confusing fast.
- You're not testing the way end users will be running *Prost*!, and hence you
  may be accidentally introducing bugs that would normally be immediately
  visible had you been running your development copy of *Prost!* outside of the
  development directory.

Caveats aside, if you wish to run *Prost!* directly from the development
directory, here's two ways:

.. code-block:: bash

    # First way: Run Prost! library module as a script (recommended)
    cd prost
    python -m prost

    # Second way: Run prost.py executable directly (less recommended)
    prost/prost.py

===================================================================
Running *Prost!* outside of the development directory (recommended)
===================================================================

Prost was designed to be installed into your PATH and run from anywhere.
Obviously during development this can become cumbersome, and so we use
setuptools to make development in this environment easier. I'll explain why in
a moment.  For now, here's the quick and dirty:

.. code-block:: bash

    # Go to Prost! codebase, and install in development mode
    cd /path/to/development/prost
    python setup.py develop
    # or
    python setup.py develop --user

    # Go to (different) directory from which you will run Prost!
    cd /some/different/directory

    # Run Prost!
    prost

    # Make changes to Prost, and re-run Prost to test changes.
    cd /path/to/development/prost
    edit prost/alignment.py ... (for example)
    cd /some/different/directory
    prost
    
    # Repeat the above block

    # All done with development for today, "uninstall" prost:
    cd /path/to/development/prost
    python setup.py develop --uninstall
    # or 
    python setup.py develop --uninstall --user

**The Why**

Instead of the details of why we do this myself, I'll instead quote this 
`great stackoverflow answer <http://sphinx.pocoo.org>`_:

    *python setup.py install* is used to install (typically third party) packages
    that you're not going to be developing/editing/debugging yourself.

    For your own stuff, you want to get your package installed and then be able to
    frequently edit your code and not have to re-install your package—this is
    exactly what python setup.py develop does: installs the package (typically just
    a source folder) in a way that allows you to conveniently edit your code after
    its installed to the (virtual) environment and have the changes take effect
    immediately.

**The Caveat**

The caveat is that running *python setup.py develop* will/may alter your
*easy-install.pth* file.  This is a system-specific file in a system-specific
path that is not always easy to find (personal experience). If you don't
realize that this file is being editing, and you have multiple development
directories, oh boy, what a **pernicious** bug to unravel. From the 
`setuptools docs
<https://pythonhosted.org/setuptools/setuptools.html#development-mode>`_:

    To do this, use the setup.py develop command. It works very similarly to
    setup.py install or the EasyInstall tool, except that it doesn’t actually
    install anything. Instead, it creates a special .egg-link file in the
    deployment directory, that links to your project’s source code. And, if your
    deployment directory is Python’s site-packages directory, it will also update
    the easy-install.pth file to include your project’s source code, thereby making
    it available on sys.path for all programs using that Python installation.

In any case, just be aware of this.  It is SUPER easy to revert your
*easy-install.pth* to the state it was before you ran *python setup.py
develop*.  All you do is run:

.. code-block:: bash

    python setup.py develop --uninstall
    # or
    python setup.py develop --uninstall --user

==============================================
Running *Prost!* from a specific tag or commit
==============================================

For performing test-data-runs (i.e. for running *Prost!* on a full dataset
where the goal is producing analyzable results, as opposed to developing and
then running *Prost!* on small test datasets), you may not want to run *Prost!*
from a branch, because that branch might change on you.  Instead you may wish
to run *Prost!* directly from a specific commit or tag, which is guaranteed not
to change.  To do so:

.. code-block:: bash
    
    # Checkout a specific tag
    git checkout tag_name

    # Checkout a specific commit
    git checkout cc92245


And then simply note the tag_name or commit SHA-1 in your lab notebook for that
particular *Prost!* run.

Note that checking out a specific tag or commit will result in the following
omimous warning:

.. code-block:: bash

    ∴ git checkout cc92245
    Note: checking out 'cc92245'.

    You are in 'detached HEAD' state. You can look around, make experimental
    changes and commit them, and you can discard any commits you make in this
    state without impacting any branches by performing another checkout.

    If you want to create a new branch to retain commits you create, you may
    do so (now or later) by using -b with the checkout command again. Example:

      git checkout -b new_branch_name

    HEAD is now at cc92245... added development environment documentation

That's ok.  As long as you only plan on running *Prost!* and not making any
changes, you have nothing to worry about.  If you do wish to make code changes,
then just follow the directions above.

When you're finished, just do:

.. code-block:: bash

    # Switch back to a working branch
    git checkout branch_name

    # Get recent changes
    git pull

===================================
Working with *Prost!* Documentation
===================================

*Prost!* documentation is quickly accessible locally under the *doc/*
directory.

To build the docs:
------------------

.. code-block:: bash

   cd prost/doc
   make html

   # Alternatively
   cd prost
   ./build.sh

To read the docs:
-----------------

Currently, the docs aren't perfectly intertwined.  Eventually, we'll have them
all linked to one another.  For now, you can access the individual pages.  Just
look for the \*.html pages.  For example, to view this document locally on OS X:

.. code-block:: bash

   cd prost/doc
   open _build/html/development.html

However, see also :ref:`sphinx_autobuild` below for an easier method.

To edit the docs:
-----------------

The docs are written in Sphinx's reStructuredText.  See
http://sphinx-doc.org/rest.html for a nice tutorial.  The docs try to follow
the `Google Python Style Guide
<http://google-styleguide.googlecode.com/svn/trunk/pyguide.html>`_, and the
parsing of this style is provided by the `Napoleon Sphinx plugin
<http://sphinx-doc.org/latest/ext/napoleon.html>`_.  This page of `Example
Google Style Python Docstrings
<http://sphinx-doc.org/latest/ext/example_google.html#example-google>`_ is
particular helpful.

Simply find the \*.rst document under the *doc/* directory that you wish to
edit, edit it, rebuild the docs, and view the result locally in your browser.

For a new document, I usually start by copying an existing \*.rst document and
modify it.

.. _sphinx_autobuild:

Using sphinx-autobuild to develop the docs:
-------------------------------------------

`Sphinx Autobuild <https://github.com/GaretJax/sphinx-autobuild>`_ is a nice
tool that runs in the background, automatically builds the docs when you save a
change, and serves the docs on a local server to your web browser at
http://127.0.0.1:8000. 

To install:

.. code-block:: bash
   
   pip install sphinx-autobuild

To use:

.. code-block:: bash
   
   # Easy
   cd docs
   make livehtml
   
   # If 'Easy' doesnt' work:
   cd /path/to/toplevel/repo
   sphinx-autobuild docs docs/_build/html/

To view, point your browser to http://127.0.0.1:8000.
