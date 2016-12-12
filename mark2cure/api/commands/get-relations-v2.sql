SELECT  `relationship_document`.`id` as `id`,
        `relationship_document`.`document_id`,
        `relationship_document`.`title` as `title`,

        COUNT(`relation_relation`.`id`) as `relation_id`,
        IF(
          SUM(`document_view`.`user_id` = {user_id})
        , true, false) as `user_completed`

FROM (
    SELECT  `document_document`.`id` as `id`,
            `document_document`.`document_id`,
            `document_document`.`title` as `title`

    FROM `document_document`

    INNER JOIN `relation_relationgroup_documents`
        ON `relation_relationgroup_documents`.`document_id` = `document_document`.`id`

    INNER JOIN `relation_relationgroup` as `group`
        ON `group`.`id` = `relation_relationgroup_documents`.`relationgroup_id`

    WHERE `group`.`enabled` = 1

) as `relationship_document`

INNER JOIN `relation_relation`
    ON `relation_relation`.`document_id` = `relationship_document`.`id`

INNER JOIN `relation_relationannotation`
    ON `relation_relationannotation`.`relation_id` = `relation_relation`.`id`

INNER JOIN `document_annotation`
    ON (`document_annotation`.`object_id` = `relation_relationannotation`.`id`
      AND `document_annotation`.`content_type_id` = 56)

INNER JOIN `document_view`
    ON `document_view`.`id` = `document_annotation`.`view_id`

GROUP BY `relationship_document`.`document_id`, `relation_relation`.`id`
