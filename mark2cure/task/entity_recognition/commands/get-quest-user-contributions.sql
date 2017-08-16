SELECT  `task_task`.`id` as `task_id`,
        `uqr`.`user_id`,
        SUM(`uqr`.`completed`) = COUNT(`task_task`.`id`) as `quest_completed`,
        COUNT(`task_task`.`id`)/COUNT(`task_task`.`id`) as `view_progress`

FROM `task_task`

INNER JOIN `task_documentquestrelationship` as `dqr`
  ON `dqr`.`task_id` = `task_task`.`id`

INNER JOIN `task_userquestrelationship` as `uqr`
  ON `uqr`.`task_id` = `task_task`.`id`

INNER JOIN `task_userquestrelationship_views` as `uqrv`
  ON `uqrv`.`userquestrelationship_id` = `uqr`.`id`
INNER JOIN `document_view`
  ON `document_view`.`id` = `uqrv`.`view_id` AND `document_view`.`task_type` = 'cr'

WHERE `task_task`.`id` = {task_id} AND `dqr`.`document_id` = {document_id}

GROUP BY `uqr`.`user_id`
