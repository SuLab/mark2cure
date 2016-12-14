SELECT  `document_document`.`id` as `doc_pk`,
        `document_document`.`document_id` as `pmid`,
        `relation_relationannotation`.`relation_id` as `relationship_id`,
        `relation_relationannotation`.`answer` as `relationship_answer`,
        `document_annotation`.`created`,
        `document_view`.`user_id`

FROM `relation_relationannotation`

INNER JOIN `document_annotation`
    ON `document_annotation`.`object_id` = `relation_relationannotation`.`id`
        AND `document_annotation`.`kind` = 'r'

INNER JOIN `document_view`
    ON `document_annotation`.`view_id` = `document_view`.`id`

INNER JOIN `document_section`
    ON `document_section`.`id` = `document_view`.`section_id`

INNER JOIN `document_document`
    ON `document_document`.`id` = `document_section`.`document_id`

WHERE `document_document`.`id` =  {document_id}
