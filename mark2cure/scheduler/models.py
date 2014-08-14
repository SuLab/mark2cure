'''
    Code adapted from
    https://groups.google.com/forum/#!msg/celery-users/CZXCh8sCK5Q/ihZgMV2HWWYJ
'''
from django.db import models
from djcelery.models import PeriodicTask, IntervalSchedule

import datetime


class TaskScheduler(models.Model):
    '''
      Used to view and inspect pending jobs (requried for "Next action will execute in X minutes")

      TODO: Learn in more detail
    '''

    periodic_task = models.ForeignKey(PeriodicTask, on_delete=models.DO_NOTHING)

    @staticmethod
    def schedule_every(task_name, period, every, args=None, kwargs=None):
        """ schedules a task by name every "every" "period". So an example call would be:
             TaskScheduler('mycustomtask', 'seconds', 30, [1,2,3])
             that would schedule your custom task to run every 30 seconds with the arguments 1,2 and 3 passed to the actual task.
        """
        permissible_periods = ['days', 'hours', 'minutes', 'seconds']
        if period not in permissible_periods:
            raise Exception('Invalid period specified')

        try:
            if int(every) < 0:
                raise Exception('Invalid "every" specified')
        except ValueError:
            raise Exception('Invalid "every" specified')

        # create the periodic task and the interval
        ptask_name = "%s_%s" % (task_name, datetime.datetime.now())  # create some name for the period task
        interval_schedules = IntervalSchedule.objects.filter(period=period, every=every)

        if interval_schedules:  # just check if interval schedules exist like that already and reuse em
            interval_schedule = interval_schedules[0]
        else:  # create a brand new interval schedule
            interval_schedule = IntervalSchedule()

            if int(every) > 0:
                interval_schedule.every = int(every)

            interval_schedule.period = period
            interval_schedule.save()
        ptask = PeriodicTask(name=ptask_name, task=task_name, interval=interval_schedule)
        if args:
            ptask.args = args
        if kwargs:
            ptask.kwargs = kwargs
        ptask.save()
        return TaskScheduler.objects.create(periodic_task=ptask)

    def stop(self):
        """ pauses the task
        """
        ptask = self.periodic_task
        ptask.enabled = False
        ptask.save()
        # Control().revoke(ptask.id, terminate=True)

    def start(self):
        ptask = self.periodic_task
        ptask.enabled = True
        ptask.save()

    def terminate(self):
        self.stop()
        ptask = self.periodic_task
        ptask.delete()

