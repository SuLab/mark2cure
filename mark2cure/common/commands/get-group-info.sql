SELECT  `group`.`id`,
        `group`.`name`,
        `group`.`stub`,
        `group`.`description`,
        COUNT(DISTINCT `task`.`id`) as `quest_count`,
        COUNT(DISTINCT `dqr`.`document_id`) as `document_count`,
        COUNT(DISTINCT `uqr`.`user_id`) as `users_engaged`,

        MIN(`uqr`.`created`) as `initial_engage_time`,
        MAX(`uqr`.`updated`) as `last_engage_time`,

        TRUNCATE(SUM(`uqr`.`completed`)/(COUNT(DISTINCT `dqr`.`document_id`)*15), 2) as `percentage_completed`

FROM `common_group` as `group`

INNER JOIN `task_task` as `task`
  ON `task`.`group_id` = `group`.`id`

INNER JOIN `task_documentquestrelationship` as  `dqr`
  ON `dqr`.`task_id` = `task`.`id`

INNER JOIN `task_userquestrelationship` as `uqr`
  ON `uqr`.`task_id` = `task`.`id`

  WHERE `group_id` = {group_pk}

GROUP BY `group`.`id`
