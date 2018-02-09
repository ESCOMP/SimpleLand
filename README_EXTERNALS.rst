Obtaining the full model code and associated scripting infrastructure
=====================================================================

#. Clone the repository. ::


#. Run the script **manage_externals/checkout_externals**. ::

      ./manage_externals/checkout_externals

   The **checkout_externals** script is a package manager that will
   populate the directory with the relevant versions of each of the
   components along with the CIME infrastructure code.

At this point you have a working version.

To see full details of how to set up a case, compile and run, see the CIME documentation at http://esmci.github.io/cime/ .

More details on checkout_externals
----------------------------------

The file **Externals.cfg** in your top-level directory tells
**checkout_externals** which tag/branch of each component should be
brought in to generate your sandbox.

**You need to rerun checkout_externals whenever Externals.cfg has
changed** (unless you have already manually updated the relevant
external(s) to have the correct branch/tag checked out). Common times
when this is needed are:

* After checking out a new branch/tag

* After merging some other branch/tag into your currently
  checked-out branch

**checkout_externals** must be run from the root of the source
tree. For example, if you cloned with::

  git clone https://github.com/escomp/ctsm.git my_ctsm_sandbox

then you must run **checkout_externals** from
``/path/to/my_ctsm_sandbox``.

To see more details of **checkout_externals**, issue ::

  ./manage_externals/checkout_externals --help

Customizing your sandbox
=============================

There are several use cases to consider when you want to customize or modify your sandbox.

Switching to a different branch or tag
-------------------------------------------

If you have already checked out a branch or tag and **HAVE NOT MADE ANY
MODIFICATIONS** it is simple to change your sandbox. Say that you
checked out ctsm1.0.0 but really wanted to have ctsm1.1.0;
you would simply do the following::

  git checkout ctsm1.1.0
  ./manage_externals/checkout_externals

You should **not** use this method if you have made any source code
changes, or if you have any ongoing cases that were created from
this sandbox. In these cases, it is often easiest to do a second **git
clone**.

Pointing to a different version of a component
----------------------------------------------

Each entry in **Externals.cfg** has the following form (we use CIME as an
example below)::

  [cime]
  local_path = cime
  protocol = git
  repo_url = https://github.com/CESM-Development/cime
  tag = cime5.4.0-alpha.20
  required = True

Each entry specifies either a tag or a branch. To point to a new tag:

#. Modify the relevant entry/entries in **Externals.cfg** (e.g., changing
   ``cime5.4.0-alpha.20`` to ``cime5.4.0-alpha.21`` above)

#. Checkout the new component(s)::

     ./manage_externals/checkout_externals

Keep in mind that changing individual components from a tag may result
in an invalid model (won't compile, won't run, not scientifically
meaningful) and is unsupported.

Committing your change to Externals.cfg
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After making this change, it's a good idea to commit the change in your
local git repository. First create a branch in your local
repository, then commit it. (Unlike with subversion, branches are stored
locally unless you explicitly push them up to github. Feel free to
create whatever local branches you'd like.) For example::

  git checkout -b my_ctsm_branch
  git add Externals.cfg
  git commit -m "Update CIME to cime5.4.0-alpha.20"

