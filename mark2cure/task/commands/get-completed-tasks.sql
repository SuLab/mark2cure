SELECT DISTINCT `document`.`id`

FROM `document_document` as `document`

LEFT JOIN `document_section` as `section`
  ON `section`.`document_id` = `document`.`id`

INNER JOIN `document_view` as `view`
  ON `view`.`section_id` = `section`.`id`
    AND `view`.`completed` = 1
    AND `view`.`user_id` = {user_id}
    AND `view`.`task_type` = '{task_type}'
