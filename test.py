import numpy as np
import MDAnalysis as mda
from MDAnalysis.transformations import fit_translation, fit_rot_trans
import MDAnalysisData
import MDAnalysis.analysis.rms
from MDAnalysis.analysis.rms import rmsd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

dataset = MDAnalysisData.membrane_peptide.fetch_membrane_peptide()

u = mda.Universe(dataset.topology, dataset.trajectory)

# compare membrane peptide movement in default trajectory
# vs. examples of fit_translation and fit_rot_trans to
# restrict movement relative to i.e., first frame
# using the same universe for both iteration and reference
# holding might be useful to expose weaknesses?

ag = u.select_atoms("protein")
ref = u.select_atoms("protein")

#----
# without any transform
native_cog = []

for ts in u.trajectory:
    print("native frame:", ts.frame)
    native_cog.append(ag.centroid())

u.trajectory[0]

#----
# remove translation in xy relative to self
# in frame 1 of same universe
transform = fit_translation(ag, ref, plane="xy",
                            weights=None)
u.trajectory.add_transformations(transform)
fit_trans_cog = []
fit_trans_rmsd = []
for ts in u.trajectory:
    print("fit_trans frame:", ts.frame)
    fit_trans_cog.append(ag.centroid())
    fit_trans_rmsd.append(rmsd(ag.positions, ref.positions))

u.trajectory[0]

#----
# same thing with fit_rot_trans
u = mda.Universe(dataset.topology, dataset.trajectory)
ag = u.select_atoms("protein")
ref_static = mda.Universe(dataset.topology, dataset.trajectory).select_atoms("protein")
transform = fit_rot_trans(ag, ref_static, plane=None,
                            weights=None)
u.trajectory.add_transformations(transform)
fit_trans_rot_cog = []
fit_trans_rot_rmsd = []
for ts in u.trajectory:
    fit_trans_rot_cog.append(ag.centroid())
    fit_trans_rot_rmsd.append(rmsd(ag.positions, ref.positions))

#----
# remove translation in xy relative to self
# in "static" second universe
u = mda.Universe(dataset.topology, dataset.trajectory)
ag = u.select_atoms("protein")
ref_static = mda.Universe(dataset.topology, dataset.trajectory).select_atoms("protein")
transform = fit_translation(ag, ref_static, plane="xy",
                            weights=None)
u.trajectory.add_transformations(transform)
static_fit_trans_cog = []
static_fit_trans_rmsd = []
for ts in u.trajectory:
    print("static fit_trans frame:", ts.frame)
    static_fit_trans_cog.append(ag.centroid())
    static_fit_trans_rmsd.append(rmsd(ag.positions, ref_static.positions))

#----
# same thing with with fit_rot_trans
u = mda.Universe(dataset.topology, dataset.trajectory)
ag = u.select_atoms("protein")
ref_static = mda.Universe(dataset.topology, dataset.trajectory).select_atoms("protein")
transform = fit_rot_trans(ag, ref_static, plane=None,
                            weights=None)
u.trajectory.add_transformations(transform)
static_fit_trans_rot_cog = []
static_fit_trans_rot_rmsd = []
for ts in u.trajectory:
    print("static fit_trans frame:", ts.frame)
    static_fit_trans_rot_cog.append(ag.centroid())
    static_fit_trans_rot_rmsd.append(rmsd(ag.positions, ref_static.positions))

#----
# remove translation in xy relative to self
# in "static" second universe, using *weights="mass"*
u = mda.Universe(dataset.topology, dataset.trajectory)
ag = u.select_atoms("protein")
print("ag.masses:", ag.masses)
ref_static = mda.Universe(dataset.topology, dataset.trajectory).select_atoms("protein")
transform = fit_translation(ag, ref_static, plane="xy",
                            weights="mass")
u.trajectory.add_transformations(transform)
static_weights_fit_trans_cog = []
static_weights_fit_trans_com = []
static_weights_fit_trans_rmsd = []
for ts in u.trajectory:
    print("weights fit_trans frame:", ts.frame)
    static_weights_fit_trans_cog.append(ag.centroid())
    static_weights_fit_trans_com.append(ag.center(weights=ag.masses))
    static_weights_fit_trans_rmsd.append(rmsd(ag.positions, ref_static.positions))
    print("ag.centroid():", ag.centroid())
    print("ag COM:", ag.center(weights=ag.masses))
    print("ref.centroid():", ref.centroid())

#----
# same thing with fit_rot_trans
u = mda.Universe(dataset.topology, dataset.trajectory)
ag = u.select_atoms("protein")
ref_static = mda.Universe(dataset.topology, dataset.trajectory).select_atoms("protein")
transform = fit_rot_trans(ag, ref_static, plane=None,
                            weights="mass")
u.trajectory.add_transformations(transform)
static_weights_fit_trans_rot_cog = []
static_weights_fit_trans_rot_com = []
static_weights_fit_trans_rot_rmsd = []
for ts in u.trajectory:
    static_weights_fit_trans_rot_cog.append(ag.centroid())
    static_weights_fit_trans_rot_com.append(ag.center(weights=ag.masses))
    static_weights_fit_trans_rot_rmsd.append(rmsd(ag.positions, ref_static.positions))

#----
# plot results
fig = plt.figure()
ax = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

for axis, data, rot_trans_data, title in zip([ax, ax2, ax3, ax4],
                             [native_cog, fit_trans_cog, static_fit_trans_cog,
                              static_weights_fit_trans_cog],
                             [native_cog, fit_trans_rot_cog, static_fit_trans_rot_cog,
                              static_weights_fit_trans_rot_cog],
                             ['native', 'fit_trans same traj',
                              'fit_trans second univ', 'with weights="mass"']):
    axis.set_title(title)
    coords = np.array(data)
    coords_rot_trans = np.array(rot_trans_data)
    axis.scatter(coords[...,0], coords[...,1], s=2, label='fit_translation')
    axis.scatter(coords_rot_trans[...,0],
                 coords_rot_trans[...,1], s=2,
                 color='brown',
                 alpha=0.2, label='fit_rot_trans')
    axis.set_xlabel('x')
    axis.set_ylabel('y')
    axis.set_xlim([0, 60])
    axis.set_ylim([0, 60])
    axis.set_xticks([0, 20, 40, 60])
    axis.set_yticks([0, 20, 40, 60])
    axis.set_aspect('equal')
    axis.legend()

fig.subplots_adjust(hspace=0.4)
fig.set_size_inches(6, 6)
fig.savefig('pr_1991.png', dpi=300)

# add another fig for COM vs centroid
# in tracking fit_trans performance
# when weights=mass is used
fig2 = plt.figure()
ax = fig2.add_subplot(111)
ax.set_title('weights="mass" AND COM tracking')
coords = np.array(static_weights_fit_trans_com)
ax.scatter(coords[...,0], coords[...,1], s=2)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_xlim([0, 60])
ax.set_ylim([0, 60])
ax.set_xticks([0, 20, 40, 60])
ax.set_yticks([0, 20, 40, 60])
ax.set_aspect('equal')
fig2.set_size_inches(6, 6)
fig2.savefig('pr_1991_fig2.png', dpi=300)
