SELECT `quest_documents`.`document_id`,
        IF(MAX(`quest_documents`.`quest_completed`) IS NULL, false, MAX(`quest_documents`.`quest_completed`)) as `quest_completed`, /* We choose the highest just incase any of them report that it was complete */

        COUNT(`document_view`.`id`) as `view_count`,
        IF(SUM(`document_view`.`completed`) = 2, true, false) as `document_view_completed`,
        IF(SUM(`document_view`.`opponent_id`) IS NULL, false, true) as `had_opponent`,
        (
          SELECT `document_pubtator`.`content`
          FROM `document_pubtator`
          WHERE `document_pubtator`.`document_id` = `quest_documents`.`document_id` AND `document_pubtator`.`kind` = 'DNorm'
          LIMIT 1
        ) as `disease`,
        (
          SELECT `document_pubtator`.`content`
          FROM `document_pubtator`
          WHERE `document_pubtator`.`document_id` = `quest_documents`.`document_id` AND `document_pubtator`.`kind` = 'GNormPlus'
          LIMIT 1
        ) as `gene`,
        (
          SELECT `document_pubtator`.`content`
          FROM `document_pubtator`
          WHERE `document_pubtator`.`document_id` = `quest_documents`.`document_id` AND `document_pubtator`.`kind` = 'tmChem'
          LIMIT 1
        ) as `drug`

FROM (
  SELECT DISTINCT `task_documentquestrelationship`.`document_id`,
         `task_documentquestrelationship`.`task_id`,
         `task_userquestrelationship`.`user_id` as `user_id`,
         IF(`task_userquestrelationship`.`id` IS NULL, false, true) as `user_quest_exists`,
         `task_userquestrelationship`.`id` as `user_quest_id`,
         `task_userquestrelationship`.`completed` as `quest_completed`

  FROM `task_documentquestrelationship`

  LEFT JOIN `task_userquestrelationship`
      ON `task_userquestrelationship`.`task_id` = `task_documentquestrelationship`.`task_id`
      AND `task_userquestrelationship`.`user_id` = {user_id}

  WHERE `task_documentquestrelationship`.`task_id` = {task_id}
  ORDER BY `task_documentquestrelationship`.`document_id` ASC
) as `quest_documents`

LEFT JOIN `task_userquestrelationship_views`
  ON `task_userquestrelationship_views`.`userquestrelationship_id` = `quest_documents`.`user_quest_id`

LEFT JOIN `document_view`
  ON `document_view`.`id` = `task_userquestrelationship_views`.`view_id`
  AND `document_view`.`task_type` = 'ner'

LEFT JOIN `document_section`
  ON `document_section`.`id` = `document_view`.`section_id`
  AND `document_section`.`document_id` = `quest_documents`.`document_id`

WHERE `quest_documents`.`user_quest_exists` = false
      OR (`quest_documents`.`user_quest_exists` = true AND `document_section`.`id` IS NOT NULL)

GROUP BY `quest_documents`.`document_id`
