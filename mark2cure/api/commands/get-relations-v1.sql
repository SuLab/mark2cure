SELECT  ANY_VALUE(`document_document`.`id`) as `id`,
        `document_document`.`document_id`,

        ANY_VALUE(`document_document`.`title`) as `title`,

        ANY_VALUE((
            SELECT COUNT(*)
            FROM `relation_relation`
            WHERE `relation_relation`.`document_id` = `document_document`.`id`
        )) as `relation_units`,

        COUNT(`view`.`id`) as `completions`,

        IF(
          SUM(`view`.`user_id` = {user_id})
        , true, false) as `user_completed`

FROM `document_document`

/* Link up the document group information for filtering purposes */
INNER JOIN `relation_relationgroup_documents`
    ON `relation_relationgroup_documents`.`document_id` = `document_document`.`id`
INNER JOIN `relation_relationgroup` as `group`
    ON `group`.`id` = `relation_relationgroup_documents`.`relationgroup_id`

/* Link up the completed relationship views for each document
   A View is considered complete when all relationships have been submitted (requires 100%) */
INNER JOIN `document_section`
    ON `document_section`.`document_id` = `document_document`.`id`

INNER JOIN `document_view` as `view`
    ON (`view`.`section_id` = `document_section`.`id`
      AND `view`.`task_type` = 'ri'
      AND `view`.`completed` = 1)

WHERE `group`.`enabled` = 1

GROUP BY `document_document`.`document_id`

/* Filter what we want to show on the dashboard */
HAVING  `relation_units` <= 20
    AND `completions` < {completions}
    AND `user_completed` = false

ORDER BY `completions` DESC

LIMIT 20
