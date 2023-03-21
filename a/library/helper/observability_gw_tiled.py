from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astroplan import (FixedTarget, Observer, AltitudeConstraint,
                       AtNightConstraint, MoonSeparationConstraint)
from astroplan.utils import time_grid_from_range
import astropy.units as u
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.time import TimeDelta
from astropy.coordinates import Angle, SkyCoord
import copy
import numpy as np
import operator
import astropy.units as u
from astropy.time import Time
from collections.abc import Sequence
import warnings
import pytz
class merged_def2():
    def _has_twin(ax):
        """
        Solution for detecting twin axes built on `ax`. Courtesy of
        Jake Vanderplas http://stackoverflow.com/a/36209590/1340208
        """
        for other_ax in ax.figure.axes:
            if other_ax is ax:
                continue
            if other_ax.bbox.bounds == ax.bbox.bounds:
                return True
        return False

    def plot_airmass_n_altitude(targets, observer, time, ax=None, style_kwargs=None,
                     style_sheet=None, brightness_shading=False,
                     altitude_yaxis=False, min_airmass=1.0, min_region=None,
                     max_airmass=3.0, max_region=None, use_local_tz=False):
        r"""
        Plots airmass as a function of time for a given target.

        If a `~matplotlib.axes.Axes` object already exists, an additional
        airmass plot will be "stacked" on it.  Otherwise, creates a new
        `~matplotlib.axes.Axes` object and plots airmass on top of that.

        When a scalar `~astropy.time.Time` object is passed in (e.g.,
        ``Time('2000-1-1')``), the resulting plot will use a 24-hour window
        centered on the time indicated, with airmass sampled at regular
        intervals throughout.
        However, the user can control the exact number and frequency of airmass
        calculations used by passing in a non-scalar `~astropy.time.Time`
        object. For instance, ``Time(['2000-1-1 23:00:00', '2000-1-1
        23:30:00'])`` will result in a plot with only two airmass measurements.

        For examples with plots, visit the documentation of
        :ref:`plots_time_dependent`.

        Parameters
        ----------
        targets : list of `~astroplan.FixedTarget` objects
            The celestial bodies of interest.
            If a single object is passed it will be converted to a list.

        observer : `~astroplan.Observer`
            The person, telescope, observatory, etc. doing the observing.

        time : `~astropy.time.Time`
            If scalar (e.g., ``Time('2000-1-1')``), will result in plotting target
            airmasses once an hour over a 24-hour window.
            If non-scalar (e.g., ``Time(['2000-1-1'])``, ``[Time('2000-1-1')]``,
            ``Time(['2000-1-1', '2000-1-2'])``),
            will result in plotting data at the exact times specified.

        ax : `~matplotlib.axes.Axes` or None, optional.
            The `~matplotlib.axes.Axes` object to be drawn on.
            If None, uses the current ``Axes``.

        style_kwargs : dict or None, optional.
            A dictionary of keywords passed into `~matplotlib.pyplot.plot_date`
            to set plotting styles.

        style_sheet : dict or `None` (optional)
            matplotlib style sheet to use. To see available style sheets in
            astroplan, print *astroplan.plots.available_style_sheets*. Defaults
            to the light theme.

        brightness_shading : bool
            Shade background of plot to scale roughly with sky brightness. Dark
            shading signifies times when the sun is below the horizon. Default
            is `False`.

        altitude_yaxis : bool
            Add alternative y-axis on the right side of the figure with target
            altitude. Default is `False`.

        min_airmass : float
            Lower limit of y-axis airmass range in the plot. Default is ``1.0``.

        max_airmass : float
            Upper limit of y-axis airmass range in the plot. Default is ``3.0``.

        min_region : float
            If set, defines an interval between ``min_airmass`` and ``min_region``
            that will be shaded. Default is `None`.

        max_region : float
            If set, defines an interval between ``max_airmass`` and ``max_region``
            that will be shaded. Default is `None`.

        use_local_tz : bool
            If the time is specified in a local timezone, the time will be plotted
            in that timezone.

        Returns
        -------
        ax : `~matplotlib.axes.Axes`
            An ``Axes`` object with added airmass vs. time plot.

        Notes
        -----
        y-axis is inverted and shows airmasses between 1.0 and 3.0 by default.
        If user wishes to change these, use ``ax.<set attribute>`` before drawing
        or saving plot:

        """
        # Import matplotlib, set style sheet
        if style_sheet is not None:
            _set_mpl_style_sheet(style_sheet)

        import matplotlib.pyplot as plt
        from matplotlib import dates

        # Set up plot axes and style if needed.
        if ax is None:
            ax = plt.gca()
        if style_kwargs is None:
            style_kwargs = {}
        style_kwargs = dict(style_kwargs)
        style_kwargs.setdefault('linestyle', '-')
        style_kwargs.setdefault('linewidth', 1.5)
        style_kwargs.setdefault('fmt', '-')

        if hasattr(time, 'utcoffset') and use_local_tz:
            tzoffset = time.utcoffset()
            tzname = time.tzname()
            tzinfo = time.tzinfo
        else:
            tzoffset = 0
            tzname = 'UTC'
            tzinfo = None
        # Populate time window if needed.
        # (plot against local time if that's requested)
        time_ut = Time(time)
        if time_ut.isscalar:
            time_ut = time_ut + np.linspace(-12, 12, 100)*u.hour
        elif len(time_ut) == 1:
            warnings.warn('You used a Time array of length 1.  You probably meant '
                          'to use a scalar. (Or maybe a list with length > 1?).',
                          PlotWarning)
        timetoplot = time_ut + tzoffset

        if not isinstance(targets, Sequence):
            targets = [targets]

        for target in targets:
            # Calculate airmass
            altitude = observer.altaz(time, target).alt
            zenith_fin = 90*u.deg-altitude
            airmass = observer.altaz(time_ut, target).secz
            # Mask out nonsense airmasses
            masked_airmass = np.ma.array(airmass, mask=airmass < 1)

            # Some checks & info for labels.
            try:
                target_name = target.name
            except AttributeError:
                target_name = ''

            # Plot data (against timezone-offset time)
            ax.plot_date(timetoplot.plot_date, masked_airmass, label=target_name, **style_kwargs)

        # Format the time axis
        xlo, xhi = (timetoplot[0]), (timetoplot[-1])
        ax.set_xlim([xlo.plot_date, xhi.plot_date])
        date_formatter = dates.DateFormatter('%H:%M')
        ax.xaxis.set_major_formatter(date_formatter)
        plt.setp(ax.get_xticklabels(), rotation=30, ha='right')

        # Shade background during night time
        if brightness_shading:
            start = time_ut[0]

            # Calculate and order twilights and set plotting alpha for each
            twilights = [
                (observer.sun_set_time(start, which='next'), 0.0),
                (observer.twilight_evening_civil(start, which='next'), 0.1),
                (observer.twilight_evening_nautical(start, which='next'), 0.2),
                (observer.twilight_evening_astronomical(start, which='next'), 0.3),
                (observer.twilight_morning_astronomical(start, which='next'), 0.4),
                (observer.twilight_morning_nautical(start, which='next'), 0.3),
                (observer.twilight_morning_civil(start, which='next'), 0.2),
                (observer.sun_rise_time(start, which='next'), 0.1),
            ]

            # add 'UTC' to each datetime object created above
            twilights = [(t[0].datetime.replace(tzinfo=pytz.utc), t[1])
                         for t in twilights]

            twilights.sort(key=operator.itemgetter(0))

            # add in left & right edges, so that if the airmass plot is requested
            # during the day, night is properly shaded
            left_edges = [(xlo.datetime.replace(tzinfo=tzinfo), twilights[0][1])] + twilights
            right_edges = twilights + [(xhi.datetime.replace(tzinfo=tzinfo), twilights[0][1])]

            for tw_left, tw_right in zip(left_edges, right_edges):
                left = tw_left[0]
                right = tw_right[0]
                if tzinfo is not None:
                    # convert to local time zone (which is plotted), then hack away the tzinfo
                    # so that matplotlib doesn't try to double down on the conversion
                    left = left.astimezone(tzinfo).replace(tzinfo=None)
                    right = right.astimezone(tzinfo).replace(tzinfo=None)
                ax.axvspan(left, right,
                           ymin=0, ymax=1, color='grey', alpha=tw_right[1])

        # Invert y-axis and set limits.
        y_lim = ax.get_ylim()
        if y_lim[1] > y_lim[0]:
            ax.invert_yaxis()
        ax.set_ylim([max_airmass, min_airmass])

        # Draw lo/hi limit regions, if present
        ymax, ymin = ax.get_ylim()       # should be (hi_limit, lo_limit)

        if max_region is not None:
            ax.axhspan(ymax, max_region, facecolor='#F9EB4E', alpha=0.10)
        if min_region is not None:
            ax.axhspan(min_region, ymin, facecolor='#F9EB4E', alpha=0.10)

        # Set labels.
        ax.set_ylabel("Airmass")
        ax.set_xlabel("Time from {0} [{1}]".format(min(timetoplot).datetime.date(), tzname))

        if altitude_yaxis and not merged_def2._has_twin(ax):
            altitude_ticks = np.array([90, 60, 50, 40, 30, 20])
            airmass_ticks = 1./np.cos(np.radians(90 - altitude_ticks))

            ax2 = ax.twinx()
            ax2.invert_yaxis()
            ax2.set_yticks(airmass_ticks)
            ax2.set_yticklabels(altitude_ticks)
            ax2.set_ylim(ax.get_ylim())
            ax2.set_ylabel('Altitude [degrees]')

        # Redraw figure for interactive sessions.
        ax.figure.canvas.draw()
        # Output.
        return ax, airmass, timetoplot, altitude, zenith_fin
    def doit(observatory, crossmatched_cat, zenith, moon_sep, hdul1, time_resolution, mode, outdir):
        # Specify observer at Keck Observatory:
        loc_observatory = Observer.at_site(str(observatory))

        # Use Sesame name resolver to get coordinates for Praesepe:
        targets=[]
        for i in range(0,len(crossmatched_cat)):
            coords = SkyCoord(ra=crossmatched_cat['RA'].values[i]*u.deg, 
                              dec=crossmatched_cat['DEC'].values[i]*u.deg)
            target=FixedTarget(coord=coords, name=str(crossmatched_cat['Tag'].values[i]))
            targets.append(target)


        # Define observing constraints:
        alt_const=90-zenith
        moon_sep_const=moon_sep
        constraints = [AtNightConstraint.twilight_astronomical(),
                       MoonSeparationConstraint(min=moon_sep_const * u.deg),
                       AltitudeConstraint(min=alt_const*u.deg, max=90*u.deg)] # ZA : 0-50 degrees
                        #zenith from alt
                        #zenith=90*u.deg-altitude
        # Define range of times to observe between
        start_time_object = Time(hdul1[1].header['DATE-OBS'], format='isot', scale='utc')
        dt2 = TimeDelta(86400, format='sec')
        stop_time_object = start_time_object + dt2
        start_time_object, stop_time_object

        start_time = start_time_object
        end_time = stop_time_object
        time_resolution = time_resolution * u.hour
        #time_resolution = 0.5 * u.hour

        # Create grid of times from ``start_time`` to ``end_time``
        # with resolution ``time_resolution``
        time_grid = time_grid_from_range([start_time, end_time],
                                         time_resolution=time_resolution)


        ################### 
        observability_grid = np.zeros((len(targets), len(time_grid)))
        c1=[]
        for i, target in enumerate(targets):
            # Evaluate each constraint
            c=(observability_grid[i, :]) = constraints[0](loc_observatory, target, times=time_grid)
            c1.append(c)
        c1=np.stack(c1, axis=0 )
        print('AtNightConstraint:')
        print(c1)
        c2=[]
        for i, target in enumerate(targets):
            # Evaluate each constraint
            c=(observability_grid[i, :])  = constraints[1](loc_observatory, target, times=time_grid)
            c2.append(c)
        c2=np.stack(c2, axis=0 )
        print('MoonSeparationConstraint:')
        print(c2)
        c3=[]
        for i, target in enumerate(targets):
            # Evaluate each constraint
            c=(observability_grid[i, :]) = constraints[2](loc_observatory, target, times=time_grid)
            c3.append(c)
        c3=np.stack(c3, axis=0 )
        print('AltitudeConstraint:')
        print(c3)
        #extent = [-0.5, -0.5+len(time_grid), -0.5, 2.5]
        extent = [-0.5, -0.5+len(time_grid), -0.5, -0.5+len(crossmatched_cat)]
        c_fin= c1 & c2  & c3
        print('MixedConstraints:')
        print(c_fin)
        ################### 
        ################### FIGURE
        if mode == 'diagnostic':
            fig, ax = plt.subplots(figsize=(10,5))
            ax.imshow(c_fin, extent=extent)

            ax.set_yticks(range(0, len(targets)))
            ax.set_yticklabels([c.name for c in targets[::-1]]) #inverted target list for sanity
            ax.set_xticks(range(len(time_grid)))
            ax.set_xticklabels([t.datetime.strftime("%H:%M") for t in time_grid])
            ax.set_xticks(np.arange(extent[0], extent[1]), minor=True)
            ax.set_yticks(np.arange(extent[2], extent[3]), minor=True)


            plt.xticks(color='black', fontsize=20)
            plt.yticks(color='black', fontsize=20)
            ax.grid(which='minor', color='grey', linestyle='-', linewidth=2)
            ax.tick_params(axis='x', which='minor', bottom='off')
            plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
            ax.tick_params(axis='y', which='minor', left='off')
            fig.subplots_adjust(left=0.35, right=0.9, top=0.9, bottom=0.1)
            outname = 'Observability_TimeGrid_@'+str(observatory)+'.pdf'
            fullname = os.path.join(outdir, outname)    
            plt.savefig(fullname)
        ################################################################
        ################################################################
        fig, ax = plt.subplots(figsize=(10,5))
        zenith_fin_l=[]
        airmass_l=[]
        altitude_l=[]
        for i, target in enumerate(targets):
            ax, airmass, timetoplot, altitude, zenith_fin=merged_def2.plot_airmass_n_altitude(target, loc_observatory, time_grid, brightness_shading=True, altitude_yaxis=True)
            zenith_fin_l.append(zenith_fin)
            airmass_l.append(airmass)
            altitude_l.append(altitude)
            if mode == 'diagnostic':
                outname = 'Observability_Atmosphere_@'+str(observatory)+'.pdf'
                fullname = os.path.join(outdir, outname) 
                plt.savefig(fullname)
        return ax, airmass_l, timetoplot, altitude_l, zenith_fin_l, c_fin, time_grid