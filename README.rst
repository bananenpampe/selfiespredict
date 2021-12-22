.. These are examples of badges you might want to add to your README:
   please update the URLs accordingly

    .. image:: https://api.cirrus-ci.com/github/<USER>/selfiespredict.svg?branch=main
        :alt: Built Status
        :target: https://cirrus-ci.com/github/<USER>/selfiespredict
    .. image:: https://readthedocs.org/projects/selfiespredict/badge/?version=latest
        :alt: ReadTheDocs
        :target: https://selfiespredict.readthedocs.io/en/stable/
    .. image:: https://img.shields.io/coveralls/github/<USER>/selfiespredict/main.svg
        :alt: Coveralls
        :target: https://coveralls.io/r/<USER>/selfiespredict
    .. image:: https://img.shields.io/pypi/v/selfiespredict.svg
        :alt: PyPI-Server
        :target: https://pypi.org/project/selfiespredict/
    .. image:: https://img.shields.io/conda/vn/conda-forge/selfiespredict.svg
        :alt: Conda-Forge
        :target: https://anaconda.org/conda-forge/selfiespredict
    .. image:: https://pepy.tech/badge/selfiespredict/month
        :alt: Monthly Downloads
        :target: https://pepy.tech/project/selfiespredict
    .. image:: https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter
        :alt: Twitter
        :target: https://twitter.com/selfiespredict
      .. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
          :alt: Project generated with PyScaffold
          :target: https://pyscaffold.org/

|

==============
SelfiesPredict
==============


    Reaction outcome prediction using SELFIES 


This Repository containes our work for the project 2 of the EPFL course CS433 in machinelearning.
We carried this project out as part of the ML4Science initiative, that allows students to join a research lab and work on a practical ML problem.
This problem was kindly provided by `Philippe Schwaller <https://pschwllr.github.io/#>`_ who also advised us over the course of the project.

In this work, we have retrained a transformer based model for chemical reaction prediction using string-based molecular representations.
We compaire the established SMILES representation with the recently developped `SELFIES <https://github.com/aspuru-guzik-group/selfies>` representation, that performs well in generative models.

In transformer based chemical reaction prediction, the prediction problem is threated as a NLP translation task, where the input or educts are threated as the language that is to be translated, and the reaction products are the desired translations.




Install notes
=============

* We recommend first creating a virtual environment::
     
     conda create --name selfies_project
     conda activate selfies_project


* The code can be installed by first cloning the repository and then running pip locally::

     git clone <link>
     cd <./cloned_repository>
     pip install -e .
     
* On Windows, the rdkit wheel might not work and git/setuptools might not be installed, for that we recommend to normal rdkit-install::
     
     conda activate selfies_project
     #make sure to uninstall the not-working pypi wheel
     pip uninstall rdkit-pypi
     conda install -c rdkit rdkit

* And then install into the environment::
        
        pip install -e .         
  
* It might be possible that the setup file has to be run seperately. Due to the limited time of the project, we were not able to identify why this is nescessary on google colab.::

        python setup.py install

.. _pyscaffold-notes:

Tests
====
* To run test, run in the selfies directory::

   python -m unittest


