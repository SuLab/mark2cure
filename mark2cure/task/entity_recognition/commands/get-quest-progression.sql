SELECT  DISTINCT `task_documentquestrelationship`.`document_id`,
        MAX(`task_userquestrelationship`.`completed`) as `parent_task_completed`,
        COUNT(`document_view`.`id`) as `view_count`,
        IF(SUM(`document_view`.`completed`) = 2, true, false) as `views_completed`,
        IF(SUM(`document_view`.`opponent_id`) IS NULL, false, true) as `had_opponent`

FROM `task_documentquestrelationship`

INNER JOIN `task_userquestrelationship`
  ON `task_userquestrelationship`.`task_id` = `task_documentquestrelationship`.`task_id`
  AND `task_userquestrelationship`.`user_id` = {user_id}

INNER JOIN `task_userquestrelationship_views`
  ON `task_userquestrelationship_views`.`userquestrelationship_id` = `task_userquestrelationship`.`id`

INNER JOIN `document_view`
  ON `document_view`.`id` = `task_userquestrelationship_views`.`view_id`
  AND `document_view`.`task_type` = 'cr'

INNER JOIN `document_section`
ON `document_section`.`id` = `document_view`.`section_id`
  AND `document_section`.`document_id` = `task_documentquestrelationship`.`document_id`

WHERE `task_documentquestrelationship`.`task_id` = {task_id}

GROUP BY `task_documentquestrelationship`.`document_id`

