Cookbook
========

This is not a section about French cuisine, although you will find pretty good recipes here and contribution from our users.

Preprocessing
--------------

Generating a .rvr file from the DB for a specific subassin
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following query retrieve records suitable for generating a .rvr file. For Turkey River at Garber (434514):

.. code-block:: sql

  WITH
  outlet_link AS (SELECT * FROM master_update WHERE link_id = 434514),
  waterhshed_links(link_id, parent_link_id) AS (SELECT master_update.link_id, master_update.parent_link FROM master_update, outlet_link WHERE master_update.left BETWEEN outlet_link.left AND outlet_link.right)
  SELECT link_id, coalesce(array_length(parents_link_id, 1), 0), parents_link_id FROM (
    SELECT a.link_id AS link_id, array_remove(array_agg(b.link_id), NULL) AS parents_link_id
    FROM waterhshed_links a LEFT JOIN waterhshed_links b ON (b.parent_link_id = a.link_id)
    GROUP BY a.link_id) a
  ORDER BY link_id;

Postprocessing
--------------

Reading the HDF5 outputs with Python
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may have to install the ``tables`` package.

.. code-block:: sh

  pip install tables

Reading `PyTables <http://www.pytables.org>`_ is pretty straighforward:

.. code-block:: python

  import tables

  # Open the file in read-only mode
  h5file = tables.open_file("outputs.h5", "r")

  # Print some info about the file
  print(h5file)

  for row in h5file.root.outputs:

      # Read current row content
      # The name of the comuns depends on your outputs settings
      link_id = row['LinkID']
      time = row['Time']
      state_0 = row['State0']

      # Do something with the values
      print "link_id:", link_id, "time:", time, "state_0:", state_0

    # Close the file
    h5file.close()

Here is an other way using the the ``numpy`` and ``h5py`` package:

.. code-block:: python

  import numpy as np
  import h5py

  # open file and iterate over each row
  with h5py.File(h5_file_path, "r") as hdf_file:
      hdf_file_content = np.array(hdf_file.get("outputs"))
      for i in range(len(hdf_file_content)):
          # read current row content
          cur_row = hdf_file_content[i]
          # get the values of each columns
          cur_linkid = cur_row[0]
          cur_time = cur_row[1]
          cur_state0 = float(cur_row[2])
